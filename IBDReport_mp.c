//==================================================================================================
// direct IBD report algorithm (implemented in C, with multi-processes - OpenMPI)
// author: Shuo Yang, sdmorrisys@gmail.com; Computational Biology Lab in Columbia University
// version: 1.0; Apr 18, 2014
//==================================================================================================
// how to compile:
// ~/mpicc -o IBDReport_mp IBDReport_mp.c -lm
// notice: "-lm" flag is used specially for pow() math library
// how to compile in c2b2 clusters: (should "qlogin" first)
// /opt/OFED/current/mpi/gcc/openmpi-1.4.2/bin/mpicc -o IBDReport_mp IBDReport_mp.c
// how to run: (after adding the routine)
// ~/mpirun -n 16 IBDReport_mp -f test_1000G... -F 1 -t 1 -e 0.01 -m 50000 -d 0 -l 100000000
//==================================================================================================
// the basic ideas for multi-processes
// 1. we can't send the filehandle to the specific process directly using OpenMPI, but rather we
//  should split the previous nexus tree file into several chunks, and find the sub-tree according to
//  the process #
// 2. we should save all the IBD segments and the tentative segments to temporary files, other than
//  keeping all of the tentative segments (meaning the trees at the boundary) at memory so we should
// think of how to make the protocol of saving in the temporary files
// 3. something we need to pass from main process to sub-process:
//  -> the number of the chromosomes
//  -> length of a tree in string format
//  we can use the collective communication routines to do this broadcasting, but it seems that we
//  should use a blocking communication --> no need at all, we can just use a p2p communication
//==================================================================================================
// other update issues from multi-processes implementation:
// 1. why direct multi-processes in C's library other than OpenMP? (should think of this)
// 2. the memory allocation for the heap (** type) can be further optimized, but little speedup gained
//  the reason that I use double pointer is I want to call the element like in a matrix -> you are stupid Shuo
//==================================================================================================
// some problems: (for multi-processes)
// 1. in clusters, before each splitted file has been completely saved, the sub-process has begun -> potential
//  read error
// 2. the intermitent kill(9) problem --> maybe from memory mesh, so allocate more than 2G memory
//==================================================================================================
// features of this multi-processes version:
// 1. CUTOFF, EPSILON, DISCRETIZATION, END command line option. END should be designated.
// 2. number of working processes selection
// 3. efficient memory usage (specially optimized for this), in-time allocation and de-allocation
// 4. input from stdin, but at this time we MUST invoke only one process
//
//
//
//==================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <sys/stat.h>
#include "mpi.h"
#define MASTER 0

long int DISCRETIZATION = 0;  // discretization number used
int IBDSEGLEN = 33;  // xxxx-xxxx:xxxxxxxxx-xxxxxxxxx;
long int CUTOFF = 0;  // cutoff used
char * END = "100000000";  // end of the recombination breakpoint
long int TREELEN = 0;  // length of a tree in string format; set this long int for fear there may be 5000 chromosomes
double RATIO = 1.2;  // used by TREELEN
int PROCESSES = 2;  // processes # used, default value is 2; this is un-useful, because we must invoke n processes in OpenMPI
long int SAMPLE = 0;  // number of chromosomes
long int TOTAL = 0;  // this is the maximum index of the matrix
char * NEXUSTREE = NULL;  // the file name in string format
FILE * out_handle = NULL;  // used for each working process to save the IBD segments and the tentative segments
char temp_dir[100] = "";
double EPSILON = 0.01;
int FORMAT = -1;  // the format of the tree file; 0 means nexus format, and 1 means newick format
int TYPE = -1;  // the input type; 1 means input form tree file, and 2 means input from stdin
//=========================== the memory space used by IBD reporting ===========================
int * mark = NULL;  // this notification mark is set for tree files combination
        // if we have not yet got the first segment, we should set this value as 0
        // we have one mark for each pair
int seq;
char * tree_whole = NULL;
long int breakpoint;
long int * breakpoint_last = NULL;
double * tMRCA_last = NULL;  // tMRCA_last and previous tMRCA_start
long int * breakpoint_start = NULL;  // use left-bottom space and right-up space to save different values
long int * co_last = NULL;  // the whole length of last IBD segment; used for input from stdin
//==============================================================================================


// some transformation function
//==================================================================================================
int string_int(char * string)  // used for the number of processes, and the tree FORMAT
{
  int i = 0;
  int sum = (int)(string[i++])-48;
  while(string[i] != '\0')
  {
    sum = sum * 10 + ((int)(string[i++]) - 48);
  }
  return(sum);
}

// for recombination breakpoint, which can be at most 100,000,000; so we use a long int. also for the sample name
long int string_long(char * string)
{
  int i = 0;
  long int sum = (long int)(string[i++])-48;
  while(string[i] != '\0')
  {
    sum = sum * 10 + ((long int)(string[i++]) - 48);
  }
  return(sum);
}

void long_string(long int num, char * string)
{
  int i = 0;
  char temp[10];
  while(num != 0)
  {
    temp[i++] = num % 10 + 48;
    num = num / 10;
  }
  temp[i] = '\0';
  int n = strlen(temp);
  for(i=0;i<n;i++)
    string[i]=temp[n-i-1];
  string[i] = '\0';
}

// for tMRCA
double string_double(char * string)
{
  int i = 0;
  double sum = 0;
  if(string[i] != '0')
  {
    sum = ((double)(string[i++])-48);
    while(string[i] != '.')
    {
      sum = sum * 10 + ((double)(string[i++]) - 48);
    }
  }
  else
  {
    i++;
  }
  i++;
  int n = i-1;
  while(string[i] != '\0' && string[i] != 'e')
  {
    sum = sum + ((double)(string[i])-48)/pow(10, (i-n));
    i++;
  }
  if(string[i] == 'e')
  {
    i++;
    if(string[i] == '-')
    {
      i++;
      n = 0;
      while(string[i] != '\0')
        n = n * 10 + (int)(string[i++] - 48);
      sum = sum / pow(10, n);
    }
    else
    {
      i++;
      n = 0;
      while(string[i] != '\0')
        n = n * 10 + (int)(string[i++] - 48);
      sum = sum * pow(10, n);
    }
  }
  return(sum);
}

// transform the matrix index to "xx-xx" string format
void name_pair(long int name, char * pair)
{
  long int sample1;
  long int sample2;
  name += 1;
  if(name % SAMPLE == 0)
  {
    sample1 = name / SAMPLE;
    sample2 = SAMPLE;
  }
  else
  {
    sample1 = name / SAMPLE + 1;
    sample2 = name % SAMPLE;
  }
  sprintf(pair, "%ld-%ld", sample1, sample2);
}
//==================================================================================================


// construct the list structure in C, and it's related operations
//==================================================================================================
typedef struct Node_int  // to store the classified samples
{
  long int element;
  struct Node_int * next;
}Node_int;

typedef struct List_int  // list of int
{
  Node_int * start;
  Node_int * end;
  int length;
}List_int;

typedef struct Node_list  // for 2-dimensional list, say, l1 here, each element pointing to Node_int
{
  List_int * start;
  struct Node_list * next;
}Node_list;

typedef struct List_list  // list of list
{
  Node_list * start;
  Node_list * end;
  int length;
}List_list;

typedef struct Node_double  // to store the tMRCA of each group of samples
{
  double element;
  struct Node_double * next;
}Node_double;

typedef struct List_double  // list of double
{
  Node_double * start;
  Node_double * end;
  int length;
}List_double;

List_list * createList_list()
{
  List_list * lHead;
  lHead = (List_list *)malloc(sizeof(List_list));
  lHead -> start = NULL;
  lHead -> end = NULL;
  lHead -> length = 0;
  return(lHead);
}

List_int * createList_int()
{
  List_int * lHead;
  lHead = (List_int *)malloc(sizeof(List_int));
  lHead -> start = NULL;
  lHead -> end = NULL;
  lHead -> length = 0;
  return(lHead);
}

List_double * createList_double()
{
  List_double * lHead;
  lHead = (List_double *)malloc(sizeof(List_double));
  lHead -> start = NULL;
  lHead -> end = NULL;
  lHead -> length = 0;
  return(lHead);
}

// we only support appending an empty list (add length 1 of course) in this function
void appendList_list(List_list * list)
{
  Node_list * node;
  node = (Node_list *)malloc(sizeof(Node_list));
  node -> start = NULL;
  node -> next = NULL;

  List_int * list1;
  list1 = (List_int *)malloc(sizeof(List_int));
  list1 -> start = NULL;
  list1 -> end = NULL;
  list1 -> length = 0;

  node -> start = list1;

  if(list->start == NULL)
  {
    list -> start = node;
    list -> end = node;
    (list -> length)++;
  }
  else
  {
    ((list -> end) -> next) = node;
    (list -> end) = node;
    (list -> length)++;
  }
}

// we only support to append an unsigned int to the last List_int in the List_list, because we only
// need this
void appendList_int(List_list * list, long int value)
{
  Node_int * node;
  node = (Node_int *)malloc(sizeof(Node_int));
  node -> element = value;
  node -> next = NULL;


  if(((list->end)->start)->start == NULL)
  {
    ((list->end)->start)->start = node;
    ((list->end)->start)->end = node;
    ((list->end)->start)->length++;
  }
  else
  {
    (((list->end)->start)->end)->next = node;
    ((list->end)->start)->end = node;
    ((list->end)->start)->length++;
  }
}

void appendList_double(List_double * list, double value)
{
  Node_double * node;
  node = (Node_double *)malloc(sizeof(Node_double));
  node -> element = value;
  node -> next = NULL;

  if(list->start == NULL)
  {
    list -> start = node;
    list -> end = node;
    (list -> length)++;
  }
  else
  {
    ((list -> end) -> next) = node;
    (list -> end) = node;
    (list -> length)++;
  }
}

void extendList_list(List_list * list)
{
  Node_list * node_list_last = list -> start;
  int i = 1;
  if(list->length > 2)
  {
    while(i < (list->length)-1)
    {
      node_list_last = node_list_last->next;
      i++;
    }
  }

  if(list->end->start->length == 0)  // not possible for this case
  {
    free(list->end->start);
    free(list->end);
    list->end = node_list_last;
    list->length -= 1;
  }
  else
  {
    if(node_list_last->start->length == 0)
    {
      free(node_list_last->start);
      node_list_last->start = list->end->start;
      free(list->end);
      list->end = node_list_last;
      list->length -= 1;
    }
    else
    {
      node_list_last->start->end->next = list->end->start->start;
      node_list_last->start->end = list->end->start->end;
      node_list_last->start->length = node_list_last->start->length + list->end->start->length;
      free(list->end->start);
      free(list->end);
      list->end = node_list_last;
      list->length -= 1;
    }
  }
}

void extendList_double(List_double * list, double value)
{
  Node_double * node_double_last = list -> start;
  int i = 1;
  if(list->length > 2)
  {
    while(i < (list->length)-1)
    {
      node_double_last = node_double_last->next;
      i++;
    }
  }

  if(node_double_last->element == 0)  // left leaf
  {
    node_double_last->element = list->end->element + value;
  }
  // right leaf
  free(list->end);
  list->end = node_double_last;
  list->length -= 1;
}

void cleanList(List_list * l1, List_double * l2)
{
  // by now, the start and end of these two lists are all pointing to the first element in those lists
  Node_int * node_int_last = l1->start->start->start;
  Node_int * node_int = node_int_last->next;
  free(node_int_last);
  while(node_int->next != NULL)  // we have more than 2 elements in this List_int
  {
    node_int_last = node_int;
    node_int = node_int_last->next;
    free(node_int_last);
  }
  free(node_int);
  free(l1->start->start);
  free(l1->start);  // the start and the end are all pointing to the same Node_list
  free(l1);

  free(l2->start);
  free(l2);
}
//==================================================================================================


// the core "pasing based direct IBD report" algorithm
//==================================================================================================
void IBD_sender(long int name, long int begin_long, long int end_long)
{
  char segment[IBDSEGLEN];
  segment[0] = '\0';  // sample1-sample2:breakpoint1-breakpoint2
  char pair[10];  // 10 is enough for even 5000 chromosomes
  pair[0] = '\0';
  name_pair(name, pair);
  sprintf(segment, "%s:%ld-%ld;\n", pair, begin_long, end_long);
  //sprintf(segment, "%ld:%ld-%ld;\n", name, begin_long, end_long);
  fputs(segment, out_handle);
}

void cal_during_par(int tree_start, long int num)
{
  //printf("%s", tree_whole + tree_start);
  //printf("%ld\n", num);

  int count = 0;
  int i = tree_start;
  List_list * l1 = createList_list();  // store the nodes and their existing list
  List_double * l2 = createList_double();  // store the MRCA of corresponding nodes

  int j = i;
  while(tree_whole[j++] != '\0');
  tree_whole[j-4] = '\0';

  // the following used only for IBD update in the while loop
  Node_list * node_list = NULL;
  Node_int * node_int = NULL;
  Node_double * node_double = NULL;
  double value;  // used by the tMRCA update

  while(tree_whole[i] != '\0')
  {
    if(tree_whole[i] == ',' || tree_whole[i] == ' ' || tree_whole[i] == ';')
    {
      i += 1;
      continue;
    }

    if(tree_whole[i] == '(')
    {
      i += 1;
      appendList_list(l1);
      appendList_double(l2, 0);  // previous None in Python
      continue;
    }

    if(tree_whole[i] >= 48 && tree_whole[i] <= 57)
    {
      char node[5];
      int j = 0;
      while(tree_whole[i] != '.')
      {
        node[j++] = tree_whole[i++];
      }
      node[j] = '\0';
      long int node_value = string_long(node);
      while(tree_whole[i++] != ':');
      char tmrca[20];  // I give this double type number 20 effective digits
      j = 0;
      while(tree_whole[i] != ',' && tree_whole[i] != ')' && tree_whole[i] != '\0')
      {
        tmrca[j++] = tree_whole[i];
        i++;
      }
      tmrca[j] = '\0';
      double tmrca_value = string_double(tmrca);

      appendList_int(l1, node_value);  // we fix this kind of appending only to the last list in the list

      if((l2->end)->element == 0)  // two sub-nodes
      {
        l2->end->element = tmrca_value;
      }

      //------- IBD update -------
      j = 0;
      node_list = l1->start;
      node_int = (node_list->start)->start;
      node_double = l2->start;

      for(j=1; j<=l1->length; j++)
      {
        do
        {
          if(node_int == NULL)  // there may be an empty list
          {
            break;
          }
          long int node1 = node_int->element;
          if(node1 != node_value)
          {
            long int name;
            if(node1 < node_value)
            {
              name = (node1 - 1) * SAMPLE + node_value - 1;  // here the name should be an long integer
            }
            else
            {
              name = (node_value - 1) * SAMPLE + node1 - 1;
            }
            //printf("%u %u\n", node1, node_value);
            if(num == 1)  // we should use the first tree to initialize the repository
            {
              breakpoint_last[name] = breakpoint;
              tMRCA_last[name] = node_double->element;

              // for the tree boundary
              if(seq != 1)
              {
                breakpoint_start[name] = breakpoint;
                tMRCA_last[TOTAL - name] = node_double->element;
              }
            }
            else
            {
              if((node_double->element - tMRCA_last[name]) > EPSILON || (node_double->element - tMRCA_last[name]) < -EPSILON)
              {//give tolerance to errors due to data transformation
                // longer enough or not
                if(seq != 1 && mark[name] == 0)
                {
                  breakpoint_start[TOTAL - name] = breakpoint;
                  breakpoint_last[name] = breakpoint;
                  tMRCA_last[name] = node_double->element;
                  mark[name] = 1;
                }
                else
                {
                  if((breakpoint - breakpoint_last[name]) > CUTOFF)
                  {
                    IBD_sender(name, breakpoint_last[name], breakpoint);
                    breakpoint_last[name] = breakpoint;
                    tMRCA_last[name] = node_double->element;
                  }
                  else
                  {
                    breakpoint_last[name] = breakpoint;
                    tMRCA_last[name] = node_double->element;
                  }
                }
              }
            }
          }
          if(node_int->next == NULL){break;}
          else{node_int = node_int->next;}
        }
        while(1);
        if(j<l1->length)
        {
          node_list = node_list->next;
          node_int = (node_list->start)->start;
          node_double = node_double->next;
        }
      }
      continue;
    }

    if(tree_whole[i] == ')')  // there must be a ":" following the ")"
    {
      // two possibilities: left leaf and right leaf
      extendList_list(l1);
      i += 2;
      char tmrca[20];  // I give this double type number 20 effective digits
      int j = 0;
      while(tree_whole[i] != ',' && tree_whole[i] != ')' && tree_whole[i] != '\0')
      {
        tmrca[j++] = tree_whole[i];
        i++;
      }
      tmrca[j] = '\0';
      value = string_double(tmrca);
      extendList_double(l2, value);
      continue;
    }
  }
  cleanList(l1, l2);
}

void cal_during_par_newick(int tree_start, long int num, long int segment)
{
  //printf("%s", tree_whole + tree_start);
  //printf("%ld\n", num);

  // For the Newick tree file (-F 2), as the given data is the length of present segment, other than
  // the beginning coordinate, we should change our algorithm:
  // 1. we have two spaces (p1, p2) for saving the temporary data, one for last coordinate, and the other
  //  for the next coordinate
  // 2. if this is the first time, simply put "0" into p1, and the segment length into p2
  // 3. otherwise, detect tMRCA changes, if change, -> 4; else -> 5
  // 4. if p2-p1> CUTOFF, we report; then, we let p1 = p2, and p2 = p2 + segment length
  // 5. we let p2 = p2 + segment length
  // 6. when we move back to the main process, we should report (or not) the last segment

  int count = 0;
  int i = tree_start;
  List_list * l1 = createList_list();  // store the nodes and their existing list
  List_double * l2 = createList_double();  // store the MRCA of corresponding nodes

  int j = i;
  while(tree_whole[j++] != '\0');
  tree_whole[j-3] = '\0';

  // the following used only for IBD update in the while loop
  Node_list * node_list = NULL;
  Node_int * node_int = NULL;
  Node_double * node_double = NULL;
  double value;  // used by the tMRCA update

  while(tree_whole[i] != '\0')
  {
    if(tree_whole[i] == ',' || tree_whole[i] == ' ' || tree_whole[i] == ';')
    {
      i += 1;
      continue;
    }

    if(tree_whole[i] == '(')
    {
      i += 1;
      appendList_list(l1);
      appendList_double(l2, 0);  // previous None in Python
      continue;
    }

    if(tree_whole[i] >= 48 && tree_whole[i] <= 57)
    {
      char node[5];
      int j = 0;
      while(tree_whole[i] != ':')
      {
        node[j++] = tree_whole[i++];
      }
      node[j] = '\0';
      long int node_value = string_long(node) + 1;
      i++;
      char tmrca[20];  // I give this double type number 20 effective digits
      j = 0;
      while(tree_whole[i] != ',' && tree_whole[i] != ')' && tree_whole[i] != '\0')
      {
        tmrca[j++] = tree_whole[i];
        i++;
      }
      tmrca[j] = '\0';
      double tmrca_value = string_double(tmrca);

      appendList_int(l1, node_value);  // we fix this kind of appending only to the last list in the list

      if((l2->end)->element == 0)  // two sub-nodes
      {
        l2->end->element = tmrca_value;
      }

      //------- IBD update -------
      j = 0;
      node_list = l1->start;
      node_int = (node_list->start)->start;
      node_double = l2->start;

      for(j=1; j<=l1->length; j++)
      {
        do
        {
          if(node_int == NULL)  // there may be an empty list
          {
            break;
          }
          long int node1 = node_int->element;
          if(node1 != node_value)
          {
            long int name;
            if(node1 < node_value)
            {
              name = (node1 - 1) * SAMPLE + node_value - 1;  // here the name should be an long integer
            }
            else
            {
              name = (node_value - 1) * SAMPLE + node1 - 1;
            }
            //printf("%u %u\n", node1, node_value);
            if(num == 1)  // we should use the first tree to initialize the repository
            {
              co_last[name] = 0;
              co_last[TOTAL - name] = segment;
              tMRCA_last[name] = node_double->element;

              // for the tree boundary
              //if(seq != 1)
              //{
              //  breakpoint_start[name] = breakpoint;
              //  tMRCA_last[TOTAL - name] = node_double->element;
              //}
            }
            else
            {
              if((node_double->element - tMRCA_last[name]) > EPSILON || (node_double->element - tMRCA_last[name]) < -EPSILON)
              {//give tolerance to errors due to data transformation
                // longer enough or not
                //if(seq != 1 && mark[name] == 0)
                //{
                //  breakpoint_start[TOTAL - name] = breakpoint;
                //  breakpoint_last[name] = breakpoint;
                //  tMRCA_last[name] = node_double->element;
                //  mark[name] = 1;
                //}
                //else
                //{
                if((co_last[TOTAL - name] - co_last[name]) > CUTOFF)
                {
                  IBD_sender(name, co_last[name], co_last[TOTAL - name]);
                  //printf("%ld:%ld,%ld\n", name, co_last[name], co_last[TOTAL - name]);
                  co_last[name] = co_last[TOTAL - name];
                  co_last[TOTAL - name] = co_last[TOTAL - name] + segment;
                  tMRCA_last[name] = node_double->element;
                }
                else
                {
                  co_last[name] = co_last[TOTAL - name];
                  co_last[TOTAL - name] = co_last[TOTAL - name] + segment;
                  tMRCA_last[name] = node_double->element;
                  }
                //}
              }
              else
              {
                co_last[TOTAL - name] = co_last[TOTAL - name] + segment;
              }
            }
          }
          if(node_int->next == NULL){break;}
          else{node_int = node_int->next;}
        }
        while(1);
        if(j<l1->length)
        {
          node_list = node_list->next;
          node_int = (node_list->start)->start;
          node_double = node_double->next;
        }
      }
      continue;
    }

    if(tree_whole[i] == ')')  // there must be a ":" following the ")"
    {
      // two possibilities: left leaf and right leaf
      extendList_list(l1);
      i += 2;
      char tmrca[20];  // I give this double type number 20 effective digits
      int j = 0;
      while(tree_whole[i] != ',' && tree_whole[i] != ')' && tree_whole[i] != '\0')
      {
        tmrca[j++] = tree_whole[i];
        i++;
      }
      tmrca[j] = '\0';
      value = string_double(tmrca);
      extendList_double(l2, value);
      continue;
    }
  }
  cleanList(l1, l2);
}
//==================================================================================================


// some functions used by main (tree-processing) function
//==================================================================================================
// judge whether a newly read string is a nexus tree
int judge()
{
  if(tree_whole[1] == 't' && tree_whole[2] == 'r' && tree_whole[3] == 'e' && tree_whole[4] == 'e')
    {return(1);}
  else
    {return(0);}
}

int judge_old(char * tree_whole)
{
  if(tree_whole[1] == 't' && tree_whole[2] == 'r' && tree_whole[3] == 'e' && tree_whole[4] == 'e')
    {return(1);}
  else
    {return(0);}
}

// extract the breakpoint and the start point of the actual tree from the whole tree string
int get_breakpoint_treestart()
{
  int i = 0;
  while(tree_whole[i] != '\0')
  {
    if(tree_whole[i] == 'p' && tree_whole[i+1] == 'o' && tree_whole[i+2] == 's' && tree_whole[i+3] == '_')
    {
      char breakpoint_char[10];
      i = i+4;
      int j = 0;
      while(tree_whole[i] != ' ')
      {
        breakpoint_char[j++] = tree_whole[i++];
      }
      breakpoint_char[j] = '\0';
      breakpoint = string_long(breakpoint_char);
      break;
    }
    i++;
  }
  while(tree_whole[i++] != '(');
  return(i-1);
}

// only the main process use this function
void tree_preprocessing(char * filename)  // we need this filename here
{
  // get the file handle
  FILE * fp;
  fp = fopen(filename, "r");

  // get the size of the original tree file
  struct stat buff;
  stat(filename, &buff);
  double size_total = buff.st_size;
  printf("File '%s' has size of %f Bytes.\n", filename, size_total);
  printf("We will split these trees into several parts according to your number of processes.\n");

  // get the size of each tree string; and the number of chromosomes
  char temp[10];
  int i = 0;
  char ch;
  while((ch=fgetc(fp))!=EOF)
  {
    //printf("%c\n", ch);
    if(ch == '\n')
    {
      i = 0;
      continue;
    }
    temp[i++] = ch;
    if(i == 5)
    {
      temp[i] = '\0';
      //printf("%s\n", temp);
      if(!strcmp((temp + 1), "tree"))
      {
        break;
      }
      else
      {
        i = 0;
        while((ch=fgetc(fp))!='\n');
        continue;
      }
    }
  }
  int count = 5;
  long int sample = 0;
  int begin = 0;
  while((ch=fgetc(fp))!='\n')
  {
    count++;
    if(ch != '(' && !begin)
      continue;
    else
      begin = 1;

    if(ch == ',' || ch == ' ' || ch == ';' || ch == '(')
    {
      continue;
    }

    if(ch >= 48 && ch <= 57)
    {
      sample += 1;
      while(ch != ',' && ch != ')')
      {
        ch = fgetc(fp);
        count++;
      }
      continue;
    }

    if(ch == ':')  // there must be a ":" following the ")"
    {
      ch = fgetc(fp);
      count += 1;
      while(ch != ',' && ch != ')')
      {
        ch = fgetc(fp);
        count++;
      }
      continue;
    }
  }
  SAMPLE = sample;
  printf("There are %ld chromosomes in this scenario.\n", SAMPLE);
  TOTAL = sample*sample - 1;
  fclose(fp);
  TREELEN = (long int)(count * RATIO);  // is this RATIO enough?
  printf("The tree_length is %ld.\n", TREELEN);

  // divide the trees
  // ---> or we can divide the trees after we scan all the trees and get its number
  int tree_num = size_total / (double)(count);
  printf("The total tree number is about %d.\n", tree_num);
  printf("We want to open %d processes for this task.\n", PROCESSES);
  printf("Each thread will deal with about %d trees.\n", (int)(tree_num * 1.0 / PROCESSES));
  int DIVISION = (int)(tree_num * 1.0 / PROCESSES);
  
  if(PROCESSES > 1)  // otherwise we have the previous tree as the working tree
  {
    // create the temporary folder
    char * predir = getcwd(NULL, 0);
    sprintf(temp_dir, "%s%s%s%s", predir, "/", NEXUSTREE, "_temp");
    mkdir(temp_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);  // with read/write/search permissions for 
                                                            // owner and group, and with read/search 
                                                            // permissions for others

    fp = fopen(filename, "r");
    FILE * fp_temp;
    int n = strlen(temp_dir) + 10;  // '/tempxxx'
    char temp_file[n];
    temp_file[0] = '\0';
    sprintf(temp_file, "%s%s", temp_dir, "/temp");
    int num = 0;
    char tree_whole[TREELEN];  // only use here temporarily, and we will use judge_old
    for(num = 1; num <= PROCESSES; num++)
    {
      char temp_file1[n];
      temp_file1[0] = '\0';
      sprintf(temp_file1, "%s%d", temp_file, (num));
      fp_temp = fopen(temp_file1, "w");
      int tree_num = 0;
      while(!feof(fp))
      {
        fgets(tree_whole, TREELEN, fp);
        if(judge_old(tree_whole))
        {
          tree_num++;
          //printf("%d\n", tree_num);
          fputs(tree_whole, fp_temp);
          if(tree_num > DIVISION && num < PROCESSES)  // not for the last sub-file
            break;
        }
      }
      fclose(fp_temp);
      printf("Temporary tree file '%s' has been generated.\n", temp_file1);
      // here, right after we get a splitted tree file, we inform that process to start
      // but in clusters, before the file has been completely saved, the sub-process has begun
      printf("That working process (#%d) will be informed and start its works right now.\n", num);
      int dest = num;
      int tag1 = 2;
      int tag2 = 1;
      MPI_Send(&TREELEN, 1, MPI_LONG, dest, tag1, MPI_COMM_WORLD);
      MPI_Send(&SAMPLE, 1, MPI_LONG, dest, tag2, MPI_COMM_WORLD);
      printf("Sent signal to task(working process) #%d\n", dest);
    }
    fclose(fp);
  }
  else
  {
    printf("The working process (#%d) will be informed and start its works right now.\n", 1);
    int dest = 1;
    int tag1 = 2;
    int tag2 = 1;
    MPI_Send(&TREELEN, 1, MPI_LONG, dest, tag1, MPI_COMM_WORLD);
    MPI_Send(&SAMPLE, 1, MPI_LONG, dest, tag2, MPI_COMM_WORLD);
    printf("Sent signal to task(working process) #%d\n", dest);
  }
}
//==================================================================================================


// functions and variables about multi-processes
//==================================================================================================
void getpackage(int number)
{
  long int i = 0;

  if(PROCESSES > 1 && number > 1)
  {
    mark = (int *)malloc(sizeof(int) * SAMPLE*SAMPLE);
      for(i=0; i<SAMPLE*SAMPLE; i++)
      {
        mark[i] = 0;
      }
    breakpoint_start = (long int *)malloc(sizeof(long int) * SAMPLE*SAMPLE);
  }

  seq = number;
  tree_whole = (char *)malloc(sizeof(char) * TREELEN);
  breakpoint_last = (long int *)malloc(sizeof(long int) * SAMPLE*SAMPLE);
  tMRCA_last = (double *)malloc(sizeof(double) * SAMPLE*SAMPLE);
}

void freepackage()
{
  long int i = 0;
  if(PROCESSES > 1 && seq > 1)
  {
    free(mark);
    free(breakpoint_start);
  }

  free(tree_whole);
  free(breakpoint_last);
  free(tMRCA_last);
}

void save_boundary()  // according to whether or not this is the first process
{
  // the nexus tree file (the sub-tree file)
  int n = strlen(temp_dir) + 15;  // + '/boundaryxxx'
  char filename[n];
  filename[0] = '\0';
  sprintf(filename, "%s%s%d", temp_dir, "/boundary", seq);
  FILE * fp;
  // I don't test whether we can or not open the nexus tree file here!!!
  //if((fp = fopen(filename, "r")) == NULL)
  //{
  //  printf("Errors when opening the nexus trees!\n");
  //  return -1;
  //}
  fp = fopen(filename, "w");

  // save all the "breakpoint_last" and "tMRCA_last" into the temporary files
  // #xxxx:xxx,xxx.xxx;
  long int i, j;
  long int name;
  char segment[50];
  for(i = 1; i <= SAMPLE; i++)
    for(j = i+1; j <= SAMPLE; j++)
    {
      name = (i - 1) * SAMPLE + j - 1;
      segment[0] = '\0';
      sprintf(segment, "#%ld:%ld,%f;", name, breakpoint_last[name], tMRCA_last[name]);
      fputs(segment, fp);
    }
  if(seq != 1)  // save all the "tMRCA_start" (in "tMRCA_last") and breakpoint_start (if mark is 0, ignore the second -> 0)
  {
    // $xxxx:xxx,xxx,xxx.xxx;
    long int i, j;
    long int name;
    char segment[50];
    for(i = 1; i <= SAMPLE; i++)
      for(j = i+1; j <= SAMPLE; j++)
      {
        name = (i - 1) * SAMPLE + j - 1;
        segment[0] = '\0';
        if(mark[name] == 0)
        {
          sprintf(segment, "$%ld:%ld,%d,%f;", name, breakpoint_start[name], 0, tMRCA_last[TOTAL - name]);
        }
        else
        {
          sprintf(segment, "$%ld:%ld,%ld,%f;", name, breakpoint_start[name], breakpoint_start[TOTAL - name], tMRCA_last[TOTAL - name]);
        }
        fputs(segment, fp);
      }
  }
  fclose(fp);
  printf("Temporary potential IBDseg file '%s' has been generated.\n", filename);
}

void tree_processing()
{
  printf("This is tree_processing process #%d (there are %d totally), and I just started.\n", seq, PROCESSES);

  if(PROCESSES == 1)
  {
    // open the saving tree file
    char * predir = getcwd(NULL, 0);
    int n = strlen(predir) + strlen(NEXUSTREE) + 10;
    char fileresult[n];
    fileresult[0] = '\0';
    sprintf(fileresult, "%s%s%s", predir, "/result_", NEXUSTREE);
    out_handle = fopen(fileresult, "w");

    // the nexus tree file (not divided, the whole)
    n = strlen(predir) + strlen(NEXUSTREE) + 2;
    char filename[n];
    filename[0] = '\0';
    sprintf(filename, "%s%s%s", predir, "/", NEXUSTREE);
    FILE * fp;
    // I don't test whether we can or not open the nexus tree file here, because we have done before.
    //if((fp = fopen(filename, "r")) == NULL)
    //{
    //  printf("Errors when opening the nexus trees!\n");
    //  return -1;
    //}
    fp = fopen(filename, "r");


    long int tree_num = 0;  // there may be many many trees
    long int lastbreakpoint = 0;
    while(!feof(fp))
    {
      fgets(tree_whole, TREELEN, fp);
      if(judge())
      {
        tree_num++;
        int tree_start = get_breakpoint_treestart();
        // judge according to the discretization
        if(tree_num == 1)
        {
          cal_during_par(tree_start, tree_num);
          lastbreakpoint = breakpoint;
        }
        else
        {
          if(breakpoint - lastbreakpoint > DISCRETIZATION)
          {
            if(tree_num % 1000 == 0)
            {
              printf("Tree# %ld in processes# %d.\n", tree_num, seq);
            }
            cal_during_par(tree_start, tree_num);
            lastbreakpoint = breakpoint;
          }
          else
          {
            tree_num--;
          }
        }
      }
    }
    // test and add/ignore the last segment to IBD report
    long int i, j;
    for(i = 1; i <= SAMPLE; i++)
      for(j = i+1; j <= SAMPLE; j++)
      {
        long int name = (i - 1) * SAMPLE + j - 1;
        if(string_long(END) - breakpoint_last[name] > CUTOFF)
        {
          IBD_sender(name, breakpoint_last[name], string_long(END));
        }
      }
    fclose(fp);
    fclose(out_handle);
  }
  else  // we have several processes
  {
    // set the temporary folder
    char * predir = getcwd(NULL, 0);
    sprintf(temp_dir, "%s%s%s%s", predir, "/", NEXUSTREE, "_temp");

    // open the saving tree file
    int n = strlen(temp_dir) + 15;  // + '/IBDtempxxx'
    char fileresult[n];
    fileresult[0] = '\0';
    sprintf(fileresult, "%s%s%d", temp_dir, "/IBDtemp", seq);
    out_handle = fopen(fileresult, "w");

    // the nexus tree file (the sub-tree file)
    n = strlen(temp_dir) + 2 + 8;  // + '/tempxxx'
    char filename[n];
    filename[0] = '\0';
    sprintf(filename, "%s%s%d", temp_dir, "/temp", seq);
    FILE * fp;
    // I don't test whether we can or not open the nexus tree file here!!!
    //if((fp = fopen(filename, "r")) == NULL)
    //{
    //  printf("Errors when opening the nexus trees!\n");
    //  return -1;
    //}
    fp = fopen(filename, "r");

    long int tree_num = 0;
    long int lastbreakpoint = 0;
    while(!feof(fp))
    {
      fgets(tree_whole, TREELEN, fp);
      if(judge())
      {
        tree_num++;
        int tree_start = get_breakpoint_treestart();
        // judge according to the discretization
        if(tree_num == 1)
        {
          cal_during_par(tree_start, tree_num);
          lastbreakpoint = breakpoint;
        }
        else
        {
          if(breakpoint - lastbreakpoint > DISCRETIZATION)
          {
            if(tree_num % 1000 == 0)
            {
              printf("Tree# %ld in process# %d.\n", tree_num, seq);
            }
            cal_during_par(tree_start, tree_num);
            lastbreakpoint = breakpoint;
          }
          else
          {
            tree_num--;
          }
        }
      }
    }
    fclose(fp);
    fclose(out_handle);
    printf("Temporary IBDseg file '%s' has been generated.\n", fileresult);
    save_boundary();  // save all the boundary trees in a separate temporary file

    // remove the temporary tree files (only for multi-processes case, not for only one process)
    remove(filename);  // success -> 0; fail -> -1; reason in errno
    printf("Temporary tree file '%s' has been removed.\n", filename);
  }
  printf("The tree_processing process #%d ends here.\n", seq);
}


void stdin_process()  // TODO by now we assume that STDIN and NEWICK happen at the same time
{
  // step: (better in one )
  // 1. set the TREELEN, get the number of chromosomes (TOTAL) from the very first tree
  // 2. processing all the trees
  PROCESSES = 1;
  TREELEN = 50000;
  tree_whole = (char *)malloc(sizeof(char) * TREELEN);
  // read until we encounter a newick tree here; and the number of chromosomes
  char temp[15] = "";
  char c;
  int i = 0;
  while(1)
  {
    // read lines from stdin
    tree_whole[0] = '\0';
    fscanf(stdin, "%[^\n]", tree_whole);
    fscanf(stdin, "%c", &c);
    if(strlen(tree_whole) < 12)
    {
      continue;
    }
    else
    {
      for(i=0;(temp[i]=tree_whole[i])!=':';i++);
      temp[i] = '\0';
      if(!strcmp(temp, "NEWICK_TREE"))
      {
        break;
      }
    }
  }
  // get the value of SAMPLE from here; in tree_whole there is a whole tree
  long int sample = 0;
  long int j = 18;  // the beginning of the actuall Newick tree
  long int len = strlen(tree_whole);
  while(j < len)
  {
    if(tree_whole[j] == ',' || tree_whole[j] == ' ' || tree_whole[j] == ';' || tree_whole[j] == '(')
    {
      j++;
      continue;
    }

    if(tree_whole[j] >= 48 && tree_whole[j] <= 57)
    {
      sample += 1;
      while(tree_whole[j] != ',' && tree_whole[j] != ')')
      {
        j++;
      }
      continue;
    }

    if(tree_whole[j] == ')')  // there must be a ":" following the ")"; or it may finish
    {
      j = j+2;
      if(j == len)
      {
        continue;
      }
      while(tree_whole[j] != ',' && tree_whole[j] != ')')
      {
        j++;
      }
      continue;
    }
  }
  SAMPLE = sample;
  printf("There are %ld chromosomes in this scenario.\n", SAMPLE);
  TOTAL = sample*sample - 1;
  // by now I have got one tree and get the value of SAMPLE

  co_last = (long int *)malloc(sizeof(long int) * SAMPLE*SAMPLE);  // coordinate last, and use this to store two coordinates
  tMRCA_last = (double *)malloc(sizeof(double) * SAMPLE*SAMPLE);

  // tree_processing || change the cal_dur_par at the same time
  // first of all, open the filehandle

  char filename[50] = "";
  char * predir = getcwd(NULL, 0);
  sprintf(filename, "%s%s%s", predir, "/", "IBDreport_stdin.txt");
  out_handle = fopen(filename, "w");


  long int tree_num = 0;
  long int lastsegment = 0;
  long int segment = 0;
  int tree_start;
  while(1)
  {

    tree_whole[0] = '\0';
    fscanf(stdin, "%[^\n]", tree_whole);
    fscanf(stdin, "%c", &c);
    if(!strcmp(tree_whole, ""))  // is this enough? because maybe the previous program is slower
    {
      break;
    }
    if(strlen(tree_whole) >= 11)
    {
      c = tree_whole[11];
      tree_whole[11] = '\0';
      if(!strcmp(tree_whole, "NEWICK_TREE"))
      {
        tree_whole[11] = c;
        // or other action from here!
        // remember to add the DISCRETIZATION here; and print the tree# every 1000 trees


        tree_num++;
        i = 14;
        temp[0] = '\0';
        while(tree_whole[i] != ']')
        {
          temp[i-14] = tree_whole[i];
          i++;
        }
        temp[i-14] = '\0';
        segment = string_long(temp);  // I get the length of this segment, long int
        while(tree_whole[i] != '(')
          i++;
        tree_start = i;  // I get the start coordinate of this tree, int
        // judge according to the discretization
        if(tree_num == 1)
        {
          cal_during_par_newick(tree_start, tree_num, segment);
          lastsegment = segment;
        }
        else
        {
          if(lastsegment > DISCRETIZATION)
          {
            if(tree_num % 1000 == 0)
            {
              printf("Tree# %ld in processes# %d.\n", tree_num, 1);
            }
            cal_during_par_newick(tree_start, tree_num, segment);
            lastsegment = segment;
          }
          else
          {
            lastsegment += segment;
            tree_num--;
          }
        }


      }
    }
  }
  // TODO here we should judge the last segment whether or not to report
  // test and add/ignore the last segment to IBD report
  long int i1, j1;
  for(i1 = 1; i1 <= SAMPLE; i1++)
    for(j1 = i+1; j1 <= SAMPLE; j1++)
    {
      long int name = (i1 - 1) * SAMPLE + j1 - 1;
      if(co_last[TOTAL - name] - co_last[name] > CUTOFF)
      {
        IBD_sender(name, co_last[name], co_last[TOTAL - name]);
      }
    }

  fclose(out_handle);
  free(tree_whole);
  free(co_last);
  free(tMRCA_last);
  printf("Processing finished for the trees from stdin.\n");
}


// I don't think put everything in global space is good. SO I just encapsulate them into a structure
// and instance this structure in heap space. So I can use a pointer to easily find them
typedef struct space
{
  long int * breakpoint_last1;
  long int * breakpoint_last2;
  double * tMRCA_last1;
  double * tMRCA_last2;
  long int * breakpoint_start;
  double * tMRCA_start;
}space;

space * get_save_space()
{
  space * s = (space *)malloc(sizeof(space));

  long int i;
  s->breakpoint_last1 = (long int *)malloc(sizeof(long int) * SAMPLE*SAMPLE);
  s->breakpoint_last2 = (long int *)malloc(sizeof(long int) * SAMPLE*SAMPLE);
  s->tMRCA_last1 = (double *)malloc(sizeof(double) * SAMPLE*SAMPLE);
  s->tMRCA_last2 = (double *)malloc(sizeof(double) * SAMPLE*SAMPLE);
  s->breakpoint_start = (long int *)malloc(sizeof(long int) * SAMPLE*SAMPLE);  // save two points
  s->tMRCA_start = (double *)malloc(sizeof(double) * SAMPLE*SAMPLE);

  return(s);
}

void free_save_space(space * s)
{
  free(s->breakpoint_last1);
  free(s->breakpoint_last2);
  free(s->tMRCA_last1);
  free(s->tMRCA_last2);
  free(s->breakpoint_start);
  free(s->tMRCA_start);
}

// only the main process uses this function
void combine_all(space * s, int number)
{
  // end: #xxxx:xxx,xxx.xxx;
  // start: $xxxx:xxx,xxx,xxx.xxx;

  //========== read the '..._IBDtempxxx' and save all the IBD segments to the final result =========
  int n = strlen(temp_dir) + 15;  // + '/IBDtempxxx'
  char filename[n];
  filename[0] = '\0';
  sprintf(filename, "%s%s%d", temp_dir, "/IBDtemp", number);
  FILE * fp = fopen(filename, "r");

  char ch;
  char segment[IBDSEGLEN];
  long int i = 0;
  while((ch=fgetc(fp))!=EOF)
  {
    if(ch != '\n')
    {
      segment[i++] = ch;
    }
    else
    {
      segment[i++] = ch;
      segment[i] = '\0'; 
      fputs(segment, out_handle);
      i = 0;
    }
  }
  fclose(fp);
  remove(filename);  // success -> 0; fail -> -1; reason in errno
  printf("Temporary IBDseg file '%s' has been removed.\n", filename);

  //========= read the '..._boundaryxxx' and save all the segment information to memory ==========
  // + '_boundaryxxx'
  filename[0] = '\0';
  sprintf(filename, "%s%s%d", temp_dir, "/boundary", number);
  fp = fopen(filename, "r");
  // end: #xxxx:xxx,xxx.xxx;
  // start: $xxxx:xxx,xxx,xxx.xxx;
  if(number == 1)  // there are only end segments for the results from first process
  {
    char ch;
    char segment[50];
    int count = 0;
    while((ch=fgetc(fp))!=EOF)
    {
      segment[count++] = ch;
      while((ch=fgetc(fp))!=';')
        segment[count++] = ch;
      segment[count++] = ch;
      segment[count] = '\0';
      count = 0;
      // save the segment to memory
      if(segment[0] == '#')  // there are only this kind of information here, and we do nothing
      {
        // save the end segments   #xxxx:xxx,xxx.xxx;
        char name_string[10];
        long int name_long;
        char breakpoint_string[10];
        long int breakpoint_long;
        char tMRCA_string[20];
        double tMRCA_double;
        int i = 1;
        int j = 0;
        while(segment[i] != ':')
          name_string[j++] = segment[i++];
        name_string[j] = '\0';
        name_long = string_long(name_string);
        j = 0;
        i++;
        while(segment[i] != ',')
          breakpoint_string[j++] = segment[i++];
        breakpoint_string[j] = '\0';
        breakpoint_long = string_long(breakpoint_string);
        j = 0;
        i++;
        while(segment[i] != ';')
          tMRCA_string[j++] = segment[i++];
        tMRCA_string[j] = '\0';
        tMRCA_double = string_double(tMRCA_string);

        s->tMRCA_last2[name_long] = tMRCA_double;
        s->breakpoint_last2[name_long] = breakpoint_long;
      }
    }
    // do nothing in this case
  }
  else
  {
    // for the present 'last' and the last 'last', we should change the pointer during each iteration
    // exchange the pointers of the present 'last' and the last 'last'
    long int * breakpoint_temp = s->breakpoint_last1;
    s->breakpoint_last1 = s->breakpoint_last2;
    s->breakpoint_last2 = breakpoint_temp;
    double * tMRCA_temp = s->tMRCA_last1;
    s->tMRCA_last1 = s->tMRCA_last2;
    s->tMRCA_last2 = tMRCA_temp;

    // extract the segment information and save them to memory until end
    char ch;
    char segment[50];
    int count = 0;
    while((ch=fgetc(fp))!=EOF)
    {
      segment[count++] = ch;
      while((ch=fgetc(fp))!=';')
        segment[count++] = ch;
      segment[count++] = ch;
      segment[count] = '\0';
      count = 0;
      // save the segment to memory
      if(segment[0] == '#')
      {
        // DEBUG
        //printf("%s\n", "####");
        // save the end segments   #xxxx:xxx,xxx.xxx;
        char name_string[10];
        long int name_long;
        char breakpoint_string[10];
        long int breakpoint_long;
        char tMRCA_string[20];
        double tMRCA_double;
        int i = 1;
        int j = 0;
        while(segment[i] != ':')
          name_string[j++] = segment[i++];
        name_string[j] = '\0';
        name_long = string_long(name_string);
        j = 0;
        i++;
        while(segment[i] != ',')
          breakpoint_string[j++] = segment[i++];
        breakpoint_string[j] = '\0';
        breakpoint_long = string_long(breakpoint_string);
        j = 0;
        i++;
        while(segment[i] != ';')
          tMRCA_string[j++] = segment[i++];
        tMRCA_string[j] = '\0';
        tMRCA_double = string_double(tMRCA_string);

        s->tMRCA_last2[name_long] = tMRCA_double;
        s->breakpoint_last2[name_long] = breakpoint_long;
        // DEBUG
        //printf("%ld\n", *(breakpoint_last2[name_long]));
      }
      if(segment[0] == '$')
      {
        //printf("%s\n", "$$$$");
        // save the start segment   $xxxx:xxx,xxx,xxx.xxx;
        char name_string[10];
        long int name_long;
        char breakpoint_string1[10];
        long int breakpoint_long1;
        char breakpoint_string2[10];
        long int breakpoint_long2;
        char tMRCA_string[20];
        double tMRCA_double;
        int i = 1;
        int j = 0;
        while(segment[i] != ':')
          name_string[j++] = segment[i++];
        name_string[j] = '\0';
        name_long = string_long(name_string);
        j = 0;
        i++;
        while(segment[i] != ',')
          breakpoint_string1[j++] = segment[i++];
        breakpoint_string1[j] = '\0';
        breakpoint_long1 = string_long(breakpoint_string1);
        j = 0;
        i++;
        while(segment[i] != ',')
          breakpoint_string2[j++] = segment[i++];
        breakpoint_string2[j] = '\0';
        breakpoint_long2 = string_long(breakpoint_string2);
        j = 0;
        i++;
        while(segment[i] != ';')
          tMRCA_string[j++] = segment[i++];
        tMRCA_string[j] = '\0';
        tMRCA_double = string_double(tMRCA_string);

        s->tMRCA_start[name_long] = tMRCA_double;
        s->breakpoint_start[name_long] = breakpoint_long1;
        s->breakpoint_start[TOTAL - name_long] = breakpoint_long2;
      }
    }
    // judge and report IBD according to the previous logic
    long int i, j;
    for(i = 1; i <= SAMPLE; i++)
      for(j = i+1; j <= SAMPLE; j++)
      {
        long int name = (i - 1) * SAMPLE + j - 1;
        //printf("%s\n", "@@");
        //printf("%ld\n", *(breakpoint_last2[name]));
        //printf("%s\n", "@@");

        // from here check all the spaces for this exact pair
        if(s->breakpoint_start[TOTAL - name] == 0 && number != PROCESSES)
        {
          // judge whether there is a tMRCA change between these two parts
          if((s->tMRCA_last1[name] - s->tMRCA_start[name]) > EPSILON || (s->tMRCA_last1[name] - s->tMRCA_start[name]) < -EPSILON)
          {
            if(s->breakpoint_start[name] - s->breakpoint_last1[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_last1[name], s->breakpoint_start[name]);
              }
          }
          else
          {
            s->breakpoint_last2[name] = s->breakpoint_last1[name];
          }
        }
        else
        {
          if(number == PROCESSES && s->breakpoint_start[TOTAL - name] == 0)
          {
            // the last part has no tMRCA change for this pair
            if((s->tMRCA_last1[name] - s->tMRCA_start[name]) > EPSILON || (s->tMRCA_last1[name] - s->tMRCA_start[name]) < -EPSILON)
            {
              if(string_long(END) - s->breakpoint_start[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_start[name], string_long(END));
              }
              if(s->breakpoint_start[name] - s->breakpoint_last1[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_last1[name], s->breakpoint_start[name]);
              }
            }
            else
            {
              if(string_long(END) - s->breakpoint_last1[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_last1[name], string_long(END));
              }
            }
          }
          else  //*(space2->mark[name]) != 0 && s != PROCESSES-2     *(space2->mark[name]) != 0 && s == PROCESSES-2
          {// normal case
            if(number == PROCESSES)
            {
              if(string_long(END) - s->breakpoint_last2[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_last2[name], string_long(END));
              }
            }
            if((s->tMRCA_last1[name] - s->tMRCA_start[name]) > EPSILON || (s->tMRCA_last1[name] - s->tMRCA_start[name]) < -EPSILON)
            {
              if(s->breakpoint_start[TOTAL - name] - s->breakpoint_start[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_start[name], s->breakpoint_start[TOTAL - name]);
              }
              if(s->breakpoint_start[name] - s->breakpoint_last1[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_last1[name], s->breakpoint_start[name]);
              }
            }
            else
            {
              if(s->breakpoint_start[TOTAL - name] - s->breakpoint_last1[name] > CUTOFF)
              {
                IBD_sender(name, s->breakpoint_last1[name], s->breakpoint_start[TOTAL - name]);
              }
            }
          }
        }
      }
  }
  fclose(fp);
  remove(filename);  // success -> 0; fail -> -1; reason in errno
  printf("Temporary potential IBDseg file '%s' has been removed.\n", filename);
  // move on to the next temporary file from next working process
}
//==================================================================================================


//================================= the entrance to the program ====================================
void input_error()
{
  printf("The parameters you entered are not in right format. Please keep the following format:\n");
  printf("./IBDReport -f TREE_FILE -F FORMAT -t TYPE -m CUTOFF -e EPSILON -d DISCRETIZATION -l LENGTH_OF_CHROMOSOME\n");
  printf("The default values are:\n");
  printf("CUTOFF: %d\n", 0);
  printf("EPSILON: %f\n", 0.01);
  printf("DISCRETIZATION: %d\n", 0);
  printf("LENGTH_OF_CHROMOSOME: %d\n", 100000000);
  printf("TREE_FILE and FORMAT and TYPE must be specified in your command line.\n");
}

int main(int argc, char * argv[])
{
  //================================ command line interaction ================================
  // format: ./IBDReport -f TREE_FILE -F FORMAT -t TYPE -m CUTOFF -e EPSILON -d DISCRETIZATION -l LENGTH_OF_CHROMOSOME
  int i;
  if(argc % 2 == 0 || argc == 1)
  {
    input_error();
    return(-1);
  }
  else
  {
    int count = 1;
    while(count < argc)
    {
      if(argv[count][0] != '-')
      {
        input_error();
        return(-1);
      }
      else
      {
        switch(argv[count][1]){
          case 'f':{
            count++;
            i = strlen(argv[count]) + 1;
            NEXUSTREE = (char *)malloc(sizeof(char) * i);
            i = 0;
            while((NEXUSTREE[i]=argv[count][i]) != '\0')
              i++;
            break;
          }
          case 'F':{
            count++;
            FORMAT = string_int(argv[count]);
            break;
          }
          case 't':{
            count++;
            TYPE = string_int(argv[count]);
            break;
          }
          case 'm':{
            count++;
            CUTOFF = string_long(argv[count]);
            break;
          }
          case 'e':{
            count++;
            EPSILON = string_double(argv[count]);
            break;
          }
          case 'd':{
            count++;
            DISCRETIZATION = string_long(argv[count]);
            break;
          }
          case 'l':{
            count++;
            int length = strlen(argv[count]) + 1;
            END = (char *)malloc(sizeof(char) * length);
            strcpy(END, argv[count]);
            break;
          }
        }
        count++;
      }
    }
    if((TYPE == 1 && NEXUSTREE == NULL) || FORMAT == -1 || TYPE == -1)
    {
      printf("The tree format (-F) and the input type (-t) are necessary in the command line.\n");
      printf("When you use tree file as input (parameter -t as 1), you must enter the tree file name.\n");
      return(-1);
    }
  }
  //========================================================================================

  if(TYPE == 2)  // from stdin, and at this time there should be only one process
  {
    struct  timeval start;
    struct  timeval end;
    double diff;
    gettimeofday(&start, NULL);

    stdin_process();

    gettimeofday(&end, NULL);

    printf("Work done for trees from stdin!...\n");
    diff = (double)(end.tv_sec-start.tv_sec)+ (double)(end.tv_usec-start.tv_usec)/1000000;
    //printf("Time used to split the original tree file is %f seconds.\n", (double)(time_finish1 - time_start)/CLOCKS_PER_SEC);
    printf("Time used totally is %f seconds.\n", diff);
  }
  else
  {
    //================== Multi-processes (OpenMPI) begin from here ===========================
    int numtasks, taskid, rc, dest, tag1, tag2, source;
    int over = 1;  // the finishing signal from each working process
    MPI_Status status;

    /***** Initializations *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    numtasks--;  // this is the actual number of working processes
    PROCESSES = numtasks;

    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    if(taskid == 0)
    {
      printf("The control process (main process) has begun.\n");
    }
    else
    {
      printf("MPI sub-task (process #%d) has opened up...\n", taskid);
    }
    tag2 = 1;
    tag1 = 2;

    // multi-processes begin from here
    if(taskid == MASTER)  // if this is the main process
    {
      // set the clock
      //clock_t time_start, time_finish1, time_finish2;
      //time_start = clock();
      struct  timeval start;
      struct  timeval end1;
      struct  timeval end2;
      double diff;
      gettimeofday(&start, NULL);

      // the nexus tree file (not divided, the whole); only try whether we can open the tree file here
      char * predir = getcwd(NULL, 0);
      int n = strlen(predir) + strlen(NEXUSTREE) + 2;
      char filename[n];
      filename[0] = '\0';
      sprintf(filename, "%s%s%s", predir, "/", NEXUSTREE);
      FILE * fp;
      if((fp = fopen(filename, "r")) == NULL)
      {
        printf("Errors when opening the nexus trees!\n");
        // close all the env of MPI
        MPI_Abort(MPI_COMM_WORLD, rc);
        return(-1);
      }
      fclose(fp);

      // get the TREELEN and the total tree number (estimated); if we want several working processes, 
      // we will divide the whole file into several small files; and get the number of chromosomes - SAMPLE
      // in multi-processes mode, when we get a sub-tree, we should start that sub-process right now
      PROCESSES = numtasks;
      tree_preprocessing(filename);

      //time_finish1 = clock();  // time used to split the nexus tree file into several chuncks
      gettimeofday(&end1, NULL);

      // waiting the working process to terminate
      if(PROCESSES == 1)
      {
        /* Wait to receive finishing signal from task #1 */
        source = 1;
        MPI_Recv(&over, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
      }
      else
      {
        printf("Now the main process is waiting for merging all the porential IBDsegs from all processes.\n");
        // opening the filehandle and allocating memory space
        n = strlen(predir) + strlen(NEXUSTREE) + 10;
        char fileresult[n];
        fileresult[0] = '\0';
        sprintf(fileresult, "%s%s%s", predir, "/result_", NEXUSTREE);
        out_handle = fopen(fileresult, "w");
        space * s = get_save_space();

        /* Wait to receive finishing signal from each task */
        for (i=1; i<=numtasks; i++)
        {
          source = i;
          MPI_Recv(&over, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
          combine_all(s, i);
          if(i == 1)
          {
            continue;
          }
          else
          {
            printf("The main process has successfully merged the boundary of #%d and #%d.\n", i, i-1);
          }
        }
        fclose(out_handle);
        free_save_space(s);
        rmdir(temp_dir);  // the directory must be empty, but now it is empty
      }

      //time_finish2 = clock();  // time used to finish all the working processes
      gettimeofday(&end2, NULL);

      // final report
      printf("Work done for tree %s!...\n", filename);
      diff = (double)(end1.tv_sec-start.tv_sec)+ (double)(end1.tv_usec-start.tv_usec)/1000000;
      //printf("Time used to split the original tree file is %f seconds.\n", (double)(time_finish1 - time_start)/CLOCKS_PER_SEC);
      printf("Time used to split the original tree file is %f seconds.\n", diff);
      diff = (double)(end2.tv_sec-end1.tv_sec)+ (double)(end2.tv_usec-end1.tv_usec)/1000000;
      //printf("Time used to get IBD for all the working processes is %f seconds.\n", (double)(time_finish2 - time_finish1)/CLOCKS_PER_SEC);
      printf("Time used to get IBD for all the working processes is %f seconds.\n", diff);
    }
    else  // otherwise this is a working process
    {
      /* Receive my start signal from the master task */
      source = MASTER;
      MPI_Recv(&TREELEN, 1, MPI_LONG, source, tag1, MPI_COMM_WORLD, &status);
      MPI_Recv(&SAMPLE, 1, MPI_LONG, source, tag2, MPI_COMM_WORLD, &status);

      // preparing and processing the sub-tree file
      PROCESSES = numtasks;
      TOTAL = SAMPLE*SAMPLE - 1;
      getpackage(taskid);
      tree_processing();

      // Send my finishing signal back to the master task
      dest = MASTER;
      MPI_Send(&over, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);

      // free the resource
      freepackage();
    }

    MPI_Finalize();
    //============================================================================================
  }

  return(0);
}