//==================================================================================================
// direct IBD report algorithm (implemented in C, with multi-threads)
// author: Shuo Yang, sdmorrisys@gmail.com; Computational Biology Lab in Columbia University
// version: 1.0; Apr 16, 2014
//==================================================================================================
// how to compile:
// gcc -o IBDReport_mt IBDReport_mt.c -lpthread -lm
//==================================================================================================
// the basic ideas for multithreads
// 1. split the previous tree files into several parts, according to the number of threads, and the number of trees;
// 2. set up several threads for each of the splitted tree files;
// 3. set up another thread to collect and save the IBD report;
// 4. pay attention to the boundary;
// 5. set signal to terminate the collecting thread;
// 6. remove all the temporary files, and end the program.
//==================================================================================================
// other update issues:
// 1. the timing is not precise
// 2. why direct multi-threads in C's library other than OpenMP? (should think of this)
// 3. the memory allocation for the heap (** type) can be further optimized, but little speedup gained
//==================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <sys/stat.h>
#include <pthread.h>

long int DISCRETIZATION = 0;  // discretization number used
int IBDSEGLEN = 33;  // xxxx-xxxx:xxxxxxxxx-xxxxxxxxx;
long int CUTOFF = 0;  // cutoff used
char * END = "100000000";  // end of the recombination breakpoint
long int TREELEN = 0;  // length of a tree in string format
double RATIO = 1.2;  // used by TREELEN
int THREADS = 6;  // threads # used
long int SAMPLE = 0;  // number of chromosomes
long int TOTAL = 0;  // this is the maximum index of the matrix
char * NEXUSTREE = NULL;  // the file name in string format

pthread_t thread[10];  // by now we support at most 9 threads for computation, one for reporting
int FINISH_TABLE[9];  // if present thread not over, 0; else, 1; after counted, 2
pthread_mutex_t mut;  // mutex used by all the threads

// used by message passing
typedef struct IBD_seg
{
	char * segment;
	struct IBD_seg * next;
}IBD_seg;

IBD_seg * entrance = NULL;
IBD_seg * report_start = NULL;

typedef struct package  // but all the contents are in heap space not stack
{
	int ** mark;  // this notification mark is set for tree files combination
					// if we have not yet got the first segment, we should set this value as 0
					// we have one mark for every pair
	int * seq;
	//char tree_whole[TREELEN];
	char * tree_whole;
	//char breakpoint[10];
	long int * breakpoint;
	//char breakpoint_last[SAMPLE*SAMPLE][10];
	long int ** breakpoint_last;
	//double tMRCA_last[SAMPLE*SAMPLE];
	double ** tMRCA_last;  // tMRCA_last and previous tMRCA_start

	//char breakpoint_start[SAMPLE*SAMPLE][10];
	//long int ** breakpoint_start1;
	//char breakpoint_start[SAMPLE*SAMPLE][10];
	//long int ** breakpoint_start2;
	// we can use only one matrix to save them
	long int ** breakpoint_start;  // use left-bottom space and right-up space to save different values
	//double tMRCA_start[SAMPLE*SAMPLE];
	//double ** tMRCA_start;  --> use "double ** tMRCA_last" to save these information
}package;
//==================================================================================================


// some transformation function
//==================================================================================================
// but SAMPLE*SAMPLE will be a problem, so we use long int for the sample size
// we have most 5,000 samples, so (unsigned int) type is enough here
/*
unsigned string_unsigned(char * string)
{
	int i = 0;
	unsigned sum = (unsigned)(string[i++])-48;
	while(string[i] != '\0')
	{
		sum = sum * 10 + ((unsigned)(string[i++]) - 48);
	}
	return(sum);
}

void unsigned_string(unsigned num, char * string)
{
	int i = 0;
	char temp[5];
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
*/
int string_int(char * string)  // used for the number of threads
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
long string_long(char * string)
{
	int i = 0;
	long int sum = (long)(string[i++])-48;
	while(string[i] != '\0')
	{
		sum = sum * 10 + ((long)(string[i++]) - 48);
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
	double sum = ((double)(string[i++])-48);
	while(string[i] != '.')
	{
		sum = sum * 10 + ((double)(string[i++]) - 48);
	}
	i++;
	int n = i-1;
	while(string[i] != '\0')
	{
		sum = sum + ((double)(string[i])-48)/pow(10, (i-n));
		i++;
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
	IBD_seg * ibd_seg = (IBD_seg *)malloc(sizeof(IBD_seg));
	ibd_seg->segment = (char *)malloc(sizeof(char) * IBDSEGLEN);
	*(ibd_seg->segment) = '\0';  // sample1-sample2:breakpoint1-breakpoint2
	ibd_seg->next = NULL;

	char pair[10];  // 10 is enough for even 5000 chromosomes
	pair[0] = '\0';
	name_pair(name, pair);
	sprintf(ibd_seg->segment, "%s:%ld-%ld;", pair, begin_long, end_long);

	// every time we have a new message, we should use 'entrance' to concatenate this message
	// ibd segment memory allocation, but we should inform the IBD_report thread to report this segment
	pthread_mutex_lock(&mut);
	// add the present segment to the reporting list
	if(report_start == NULL && entrance == NULL)
	{
		report_start = ibd_seg;
		entrance = ibd_seg;
	}
	else
	{
		entrance->next = ibd_seg;
		entrance = ibd_seg;
	}
	pthread_mutex_unlock(&mut);
}

void cal_during_par(package * space, int tree_start, long int num)
{
	//printf("%s", tree_whole + tree_start);
	//printf("%d\n", num);

	int count = 0;
	int i = tree_start;
	List_list * l1 = createList_list();  // store the nodes and their existing list
	List_double * l2 = createList_double();  // store the MRCA of corresponding nodes

	int j = i;
	while(space->tree_whole[j++] != '\0');
	space->tree_whole[j-4] = '\0';

	// the following used only for IBD update in the while loop
	Node_list * node_list = NULL;
	Node_int * node_int = NULL;
	Node_double * node_double = NULL;
	double value;  // used by the tMRCA update

	while(space->tree_whole[i] != '\0')
	{
		if(space->tree_whole[i] == ',' || space->tree_whole[i] == ' ' || space->tree_whole[i] == ';')
		{
			i += 1;
			continue;
		}

		if(space->tree_whole[i] == '(')
		{
			i += 1;
			appendList_list(l1);
			appendList_double(l2, 0);  // previous None in Python
			continue;
		}

		if(space->tree_whole[i] >= 48 && space->tree_whole[i] <= 57)
		{
			char node[5];
			int j = 0;
			while(space->tree_whole[i] != '.')
			{
				node[j++] = space->tree_whole[i++];
			}
			node[j] = '\0';
			long int node_value = string_long(node);
			while(space->tree_whole[i++] != ':');
			char tmrca[20];  // I give this double type number 20 effective digits
			j = 0;
			while(space->tree_whole[i] != ',' && space->tree_whole[i] != ')' && space->tree_whole[i] != '\0')
			{
				tmrca[j++] = space->tree_whole[i];
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
							*(space->breakpoint_last[name]) = *(space->breakpoint);
							*(space->tMRCA_last[name]) = node_double->element;

							// for the tree boundary
							if(*(space->seq) != 1)
							{
								*(space->breakpoint_start[name]) = *(space->breakpoint);
								*(space->tMRCA_last[TOTAL - name]) = node_double->element;
							}
						}
						else
						{
							if((node_double->element - *(space->tMRCA_last[name])) > 0.01 || (node_double->element - *(space->tMRCA_last[name])) < -0.01)
							{//give tolerance to errors due to data transformation
								// longer enough or not
								if(*(space->seq) != 1 && *(space->mark[name]) == 0)  // lazy evaluation, so it's safe
								{
									*(space->breakpoint_start[TOTAL - name]) = *(space->breakpoint);
									*(space->breakpoint_last[name]) = *(space->breakpoint);
									*(space->tMRCA_last[name]) = node_double->element;
									*(space->mark[name]) = 1;
								}
								else
								{
									if((*(space->breakpoint) - *(space->breakpoint_last[name])) > CUTOFF)
									{
										IBD_sender(name, *(space->breakpoint_last[name]), *(space->breakpoint));
										*(space->breakpoint_last[name]) = *(space->breakpoint);
										*(space->tMRCA_last[name]) = node_double->element;
									}
									else
									{
										*(space->breakpoint_last[name]) = *(space->breakpoint);
										*(space->tMRCA_last[name]) = node_double->element;
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

		if(space->tree_whole[i] == ')')  // there must be a ":" following the ")"
		{
			// two possibilities: left leaf and right leaf
			extendList_list(l1);
			i += 2;
			char tmrca[20];  // I give this double type number 20 effective digits
			int j = 0;
			while(space->tree_whole[i] != ',' && space->tree_whole[i] != ')' && space->tree_whole[i] != '\0')
			{
				tmrca[j++] = space->tree_whole[i];
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


// some functions used by main function
//==================================================================================================
// judge whether a newly read string is a nexus tree
int judge(package * space)
{
	if(space->tree_whole[1] == 't' && space->tree_whole[2] == 'r' && space->tree_whole[3] == 'e' && space->tree_whole[4] == 'e')
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
int get_breakpoint_treestart(package * space)
{
	int i = 0;
	while(space->tree_whole[i] != '\0')
	{
		if(space->tree_whole[i] == 'p' && space->tree_whole[i+1] == 'o' && space->tree_whole[i+2] == 's' && space->tree_whole[i+3] == '_')
		{
			char breakpoint[10];
			i = i+4;
			int j = 0;
			while(space->tree_whole[i] != ' ')
			{
				breakpoint[j++] = space->tree_whole[i++];
			}
			breakpoint[j] = '\0';
			*(space->breakpoint) = string_long(breakpoint);
			break;
		}
		i++;
	}
	while(space->tree_whole[i++] != '(');
	return(i-1);
}

void tree_preprocessing(char * filename)
{
	// get the file handle
	FILE * fp;
	fp = fopen(filename, "r");

	// get the size of the original tree file
    struct stat buff;
    stat(filename, &buff);
    double size_total = buff.st_size;
    printf("File '%s' has size of %f Bytes.\n", filename, size_total);
    printf("We will split these trees into several parts according to your number of threads.\n");

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
    printf("We want to open %d threads for this task.\n", THREADS);
    printf("Each thread will deal with about %d trees.\n", (int)(tree_num * 1.0 / THREADS));
    int DIVISION = (int)(tree_num * 1.0 / THREADS);
    
    if(THREADS > 1)  // otherwise we have the previous tree as the working tree
    {
		fp = fopen(filename, "r");
		FILE * fp_temp;
	    int n = strlen(filename) + 8;  // '_temp1'
	    char temp_file[n];
	    i = 0;
	    while((temp_file[i]=filename[i++])!='\0');
	    strcat(temp_file, "_temp");
	    int num = 0;
	    char tree_whole[TREELEN];  // only use here temporarily, and we will use judge_old
		for(num = 0; num < THREADS; num++)
	    {
	    	temp_file[i+4] = (char)(num+49);
	    	temp_file[i+5] = '\0';
	    	fp_temp = fopen(temp_file, "w");
	    	int tree_num = 0;
			while(!feof(fp))
			{
				fgets(tree_whole, TREELEN, fp);
				if(judge_old(tree_whole))
				{
					tree_num++;
					//printf("%d\n", tree_num);
					fputs(tree_whole, fp_temp);
					if(tree_num > DIVISION && num < THREADS-1)  // not for the last sub-file
						break;
				}
				// here we judge according to the descretization, and we report IBD
			}
			fclose(fp_temp);
			printf("Temporary tree file '%s' has been generated.\n", temp_file);
	    }
	    fclose(fp);
	    // remember to remove all the temporary tree files finally
    }
}
//==================================================================================================


// functions and variables about multithreads
//==================================================================================================
void getpackage(package * p, int number)
{
	long int i = 0;

	if(THREADS > 1 && number > 1)
	{
		//int mark[SAMPLE*SAMPLE] for higher threads;
		p->mark = (int **)malloc(sizeof(int *) * SAMPLE*SAMPLE);
	    for(i=0; i<SAMPLE*SAMPLE; i++)
	    {
	    	p->mark[i] = (int *)malloc(sizeof(int));
	    	*(p->mark[i]) = 0;
	    }
	}

	p->seq = (int *)malloc(sizeof(int));
	*(p->seq) = number;
	//char tree_whole[25000];
	p->tree_whole = (char *)malloc(sizeof(char) * TREELEN);
	//char breakpoint[10];
	p->breakpoint = (long int *)malloc(sizeof(long int));
	//char breakpoint_last[SAMPLE*SAMPLE][10];
	p->breakpoint_last = (long int **)malloc(sizeof(long int *) * SAMPLE*SAMPLE);
    for(i=0; i<SAMPLE*SAMPLE; i++)
    	p->breakpoint_last[i] = (long int *)malloc(sizeof(long int));
	//double tMRCA_last[SAMPLE*SAMPLE];
	p->tMRCA_last = (double **)malloc(sizeof(double *) * SAMPLE*SAMPLE);
    for(i=0; i<SAMPLE*SAMPLE; i++)
    	p->tMRCA_last[i] = (double *)malloc(sizeof(double));

    if(number != 1)
    {
	    //char breakpoint_start[SAMPLE*SAMPLE][10];
		//p->breakpoint_start1 = (long int **)malloc(sizeof(long int *) * SAMPLE*SAMPLE);
	    //for(i=0; i<SAMPLE*SAMPLE; i++)
	    //	p->breakpoint_start1[i] = (long int *)malloc(sizeof(long int));
	    //char breakpoint_start[SAMPLE*SAMPLE][10];
		//p->breakpoint_start2 = (long int **)malloc(sizeof(long int *) * SAMPLE*SAMPLE);
	    //for(i=0; i<SAMPLE*SAMPLE; i++)
	    //	p->breakpoint_start2[i] = (long int *)malloc(sizeof(long int));
		//double tMRCA_start[SAMPLE*SAMPLE];

	    p->breakpoint_start = (long int **)malloc(sizeof(long int *) * SAMPLE*SAMPLE);
	    for(i=0; i<SAMPLE*SAMPLE; i++)
	    	p->breakpoint_start[i] = (long int *)malloc(sizeof(long int));

		//p->tMRCA_start = (double **)malloc(sizeof(double *) * SAMPLE*SAMPLE);  
	    //for(i=0; i<SAMPLE*SAMPLE; i++)
	    //	p->tMRCA_start[i] = (double *)malloc(sizeof(double));
    }
    else
    {
    	//p->breakpoint_start1 = NULL;
    	//p->breakpoint_start2 = NULL;
    	p->breakpoint_start = NULL;
    	//p->tMRCA_start = NULL;

    }
}

void freepackage(package * p)
{
	long int i = 0;
	if(*(p->seq) != 1)
    {
	    //char breakpoint_start[SAMPLE*SAMPLE][10];
	    //for(i=0; i<SAMPLE*SAMPLE; i++)
	    //	free(p->breakpoint_start1[i]);
	    //free(p->breakpoint_start1);
	    //char breakpoint_start[SAMPLE*SAMPLE][10];
	    //for(i=0; i<SAMPLE*SAMPLE; i++)
	    //	free(p->breakpoint_start2[i]);
	    //free(p->breakpoint_start2);
	    for(i=0; i<SAMPLE*SAMPLE; i++)
	    	free(p->breakpoint_start[i]);
	    free(p->breakpoint_start);


		//double tMRCA_start[SAMPLE*SAMPLE];
	    //for(i=0; i<SAMPLE*SAMPLE; i++)
	    //	free(p->tMRCA_start[i]);
	    //free(p->tMRCA_start);
    }

    if(THREADS > 1 && *(p->seq) > 1)
    {
	    //int mark[SAMPLE*SAMPLE];
	    for(i=0; i<SAMPLE*SAMPLE; i++)
	    	free(p->mark[i]);
	    free(p->mark);
    }

	free(p->seq);
	//char tree_whole[25000];
	free(p->tree_whole);
	//char breakpoint[10];
	free(p->breakpoint);
	//char breakpoint_last[SAMPLE*SAMPLE][10];
    for(i=0; i<SAMPLE*SAMPLE; i++)
    	free(p->breakpoint_last[i]);
    free(p->breakpoint_last);
	//double tMRCA_last[SAMPLE*SAMPLE];
    for(i=0; i<SAMPLE*SAMPLE; i++)
    	free(p->tMRCA_last[i]);
    free(p->tMRCA_last);
}

int IBD_report_over = 0;  // to inform (from the main thread) this thread to terminate

void * IBD_report(void * fp_result)
{
	printf("This is IBD_report thread, and I just started.\n");

	FILE * out_handle = (FILE *)fp_result;
	// should monitor the 'IBD_seg * entrance' to enter the message heap, then monitor the *next to go on
	while(report_start == NULL);
	IBD_seg * pointer = report_start;
	IBD_seg * pointer_new = NULL;  // used as an intermediate pointer to free the allocated heap memory
	fputs(pointer->segment, out_handle);
	fputs("\n", out_handle);
	while(!IBD_report_over || pointer->next != NULL)
	{
		if(pointer->next == NULL)
			continue;
		pointer_new = pointer->next;
		free(pointer->segment);
		free(pointer);
		pointer = pointer_new;
		fputs(pointer->segment, out_handle);
		fputs("\n", out_handle);
	}
	free(pointer->segment);
	free(pointer);

	// check IBD_report_over to determine whether or not we should terminate this thread, but the *next should be NULL
	printf("The IBD_report thread ends here.\n");
	pthread_exit(NULL);
}

void * tree_processing(void * p)
{
	package * space = (package *)p;
	printf("This is tree_processing thread #%d (there are %d totally), and I just started.\n", *(space->seq), THREADS);

	if(THREADS == 1)
	{
		// the nexus tree file (not divided, the whole)
		char * predir = getcwd(NULL, 0);
	    int n = strlen(predir) + strlen(NEXUSTREE) + 2;
	    char filename[n];
	    long int i = 0;
	    while((filename[i]=predir[i++])!='\0');
	    strcat(filename, "/");
	    strcat(filename, NEXUSTREE);
		FILE * fp;
		// I don't test whether we can or not open the nexus tree file here, because we have done before.
		//if((fp = fopen(filename, "r")) == NULL)
		//{
		//	printf("Errors when opening the nexus trees!\n");
		//	return -1;
		//}
		fp = fopen(filename, "r");

		long int tree_num = 0;  // there may be many many trees
		long int lastbreakpoint = 0;
		while(!feof(fp))
		{
			fgets(space->tree_whole, TREELEN, fp);
			if(judge(space))
			{
				tree_num++;
				int tree_start = get_breakpoint_treestart(space);
				// judge according to the discretization
				if(tree_num == 1)
				{
					cal_during_par(space, tree_start, tree_num);
					lastbreakpoint = *(space->breakpoint);
				}
				else
				{
					if(*(space->breakpoint) - lastbreakpoint > DISCRETIZATION)
					{
						if(tree_num % 1000 == 0)
						{
							printf("Tree# %ld in thread# %d.\n", tree_num, *(space->seq));
						}
						cal_during_par(space, tree_start, tree_num);
						lastbreakpoint = *(space->breakpoint);
					}
					else
					{
						tree_num--;
					}
				}
			}
		}

		// test and add/ignore the last segment to IBD report
		long int j;
		for(i = 1; i <= SAMPLE; i++)
			for(j = i+1; j <= SAMPLE; j++)
			{
				long int name = (i - 1) * SAMPLE + j - 1;
				if(string_long(END) - *(space->breakpoint_last[name]) > CUTOFF)
				{
					IBD_sender(name, *(space->breakpoint_last[name]), string_long(END));
				}
				continue;
			}
		fclose(fp);
	}
	else
	{
		// the nexus tree file (the sub-tree file)
		char * predir = getcwd(NULL, 0);
	    int n = strlen(predir) + strlen(NEXUSTREE) + 2 + 6;  // + '_tempx'
	    char filename[n];
	    int i = 0;
	    while((filename[i]=predir[i++])!='\0');
	    strcat(filename, "/");
	    strcat(filename, NEXUSTREE);
	    strcat(filename, "_temp");
	    char s_temp[2];
	    s_temp[0] = (char)(*(space->seq)+48);
	    s_temp[1] = '\0';
	    strcat(filename, s_temp);

		FILE * fp;
		// I don't test whether we can or not open the nexus tree file here!!!
		//if((fp = fopen(filename, "r")) == NULL)
		//{
		//	printf("Errors when opening the nexus trees!\n");
		//	return -1;
		//}
		fp = fopen(filename, "r");

		long int tree_num = 0;
		long lastbreakpoint = 0;
		while(!feof(fp))
		{
			fgets(space->tree_whole, TREELEN, fp);
			if(judge(space))
			{
				tree_num++;
				int tree_start = get_breakpoint_treestart(space);
				// judge according to the discretization
				if(tree_num == 1)
				{
					cal_during_par(space, tree_start, tree_num);
					lastbreakpoint = *(space->breakpoint);
				}
				else
				{
					if(*(space->breakpoint) - lastbreakpoint > DISCRETIZATION)
					{
						if(tree_num % 1000 == 0)
						{
							printf("Tree# %ld in thread# %d.\n", tree_num, *(space->seq));
						}
						cal_during_par(space, tree_start, tree_num);
						lastbreakpoint = *(space->breakpoint);
					}
					else
					{
						tree_num--;
					}
				}
			}
		}
		fclose(fp);
		// to inform the global env that I am done.
		FINISH_TABLE[*(space->seq)-1] = 1;
	}
	printf("The tree_processing thread #%d ends here.\n", *(space->seq));
	pthread_exit(NULL);
}

void combine(package * p1, package * p2)
{
	// here we should combine these two memory space: p[num] and p[num+1]
	package * space1 = p1;
	package * space2 = p2;
	long int i, j;
	for(i = 1; i <= SAMPLE; i++)
		for(j = i+1; j <= SAMPLE; j++)
		{
			long int name = (i - 1) * SAMPLE + j - 1;
			if((*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) > 0.01 || (*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) < -0.01)
			{
				if(*(space2->breakpoint_start[TOTAL - name]) - *(space2->breakpoint_start[name]) > CUTOFF)
				{
					IBD_sender(name, *(space2->breakpoint_start[name]), *(space2->breakpoint_start[TOTAL - name]));
				}
				if(*(space2->breakpoint_start[name]) - *(space1->breakpoint_last[name]) > CUTOFF)
				{
					IBD_sender(name, *(space1->breakpoint_last[name]), *(space2->breakpoint_start[name]));
				}
			}
			else
			{
				if(*(space2->breakpoint_start[TOTAL - name]) - *(space1->breakpoint_last[name]) > CUTOFF)
				{
					IBD_sender(name, *(space1->breakpoint_last[name]), *(space2->breakpoint_start[TOTAL - name]));
				}
			}
			continue;
		}
}

void combine_all(package p[])
{
	long int i, j;
	for(i = 1; i <= SAMPLE; i++)
		for(j = i+1; j <= SAMPLE; j++)
		{
			long int name = (i - 1) * SAMPLE + j - 1;

			// from here check all the spaces for this exact pair
			int s = 0;
			while(s<THREADS-1)
			{
				package * space1 = p+s;
				package * space2 = p+s+1;

				if(*(space2->mark[name]) == 0 && s != THREADS-2)
				{
					// judge whether there is a tMRCA change between these two parts
					if((*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) > 0.01 || (*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) < -0.01)
					{
						if(*(space2->breakpoint_start[name]) - *(space1->breakpoint_last[name]) > CUTOFF)
							{
								IBD_sender(name, *(space1->breakpoint_last[name]), *(space2->breakpoint_start[name]));
							}
					}
					else
					{
						*(space2->breakpoint_last[name]) = *(space1->breakpoint_last[name]);
					}
				}
				else
					if(s == THREADS-2 && *(space2->mark[name]) == 0)
					{
						// the last part has no tMRCA change for this pair
						if((*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) > 0.01 || (*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) < -0.01)
						{
							if(string_long(END) - *(space2->breakpoint_start[name]) > CUTOFF)
							{
								IBD_sender(name, *(space2->breakpoint_start[name]), string_long(END));
							}
							if(*(space2->breakpoint_start[name]) - *(space1->breakpoint_last[name]) > CUTOFF)
							{
								IBD_sender(name, *(space1->breakpoint_last[name]), *(space2->breakpoint_start[name]));
							}
						}
						else
						{
							if(string_long(END) - *(space1->breakpoint_last[name]) > CUTOFF)
							{
								IBD_sender(name, *(space1->breakpoint_last[name]), string_long(END));
							}
						}
					}
					else  //*(space2->mark[name]) != 0 && s != THREADS-2     *(space2->mark[name]) != 0 && s == THREADS-2
					{// normal case
						if(s == THREADS-2)
						{
							if(string_long(END) - *(space2->breakpoint_last[name]) > CUTOFF)
							{
								IBD_sender(name, *(space2->breakpoint_last[name]), string_long(END));
							}
						}
						if((*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) > 0.01 || (*(space1->tMRCA_last[name]) - *(space2->tMRCA_last[TOTAL - name])) < -0.01)
						{
							if(*(space2->breakpoint_start[TOTAL - name]) - *(space2->breakpoint_start[name]) > CUTOFF)
							{
								IBD_sender(name, *(space2->breakpoint_start[name]), *(space2->breakpoint_start[TOTAL - name]));
							}
							if(*(space2->breakpoint_start[name]) - *(space1->breakpoint_last[name]) > CUTOFF)
							{
								IBD_sender(name, *(space1->breakpoint_last[name]), *(space2->breakpoint_start[name]));
							}
						}
						else
						{
							if(*(space2->breakpoint_start[TOTAL - name]) - *(space1->breakpoint_last[name]) > CUTOFF)
							{
								IBD_sender(name, *(space1->breakpoint_last[name]), *(space2->breakpoint_start[TOTAL - name]));
							}
						}
					}
				s++;
			}
		}
}
//==================================================================================================


//================================= the entrance to the program ====================================
void input_error()
{
	printf("The parameters you entered are not in right format. Please keep the following format:\n");
	printf("./IBDReport -f TREE_FILE -m CUTOFF -d DISCRETIZATION -l LENGTH_OF_CHROMOSOME -t NUMBER_OF_THREADS\n");
	printf("The default values are:\n");
	printf("CUTOFF: %d\n", 0);
	printf("DISCRETIZATION: %d\n", 0);
	printf("LENGTH_OF_CHROMOSOME: %d\n", 100000000);
	printf("NUMBER_OF_THREADS: %d\n", 6);
}

int main(int argc, const char * argv[])
{
	// here we should add a command line option for user to select parameters of the running program
	// format: ./IBDReport -f TREE_FILE -m CUTOFF -d DISCRETIZATION -l LENGTH_OF_CHROMOSOME -t NUMBER_OF_THREADS

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
						while((NEXUSTREE[i]=argv[count][i++]) != '\0');
						break;
					}
					case 'm':{
						count++;
						CUTOFF = string_long((char *)argv[count]);
						break;
					}
					case 'd':{
						count++;
						DISCRETIZATION = string_long((char *)argv[count]);
						break;
					}
					case 'l':{
						count++;
						int length = strlen(argv[count]) + 1;
						END = (char *)malloc(sizeof(char) * length);
						strcpy(END, argv[count]);
						break;
					}
					case 't':{
						count++;
						THREADS = string_int((char *)argv[count]);
						break;
					}
				}
				count++;
			}
		}
		if(NEXUSTREE == NULL)
		{
			printf("There must be a tree file name existing in the input parameters.\n");
			return(-1);
		}
	}

	// set the clock
	//clock_t time_start, time_finish;
	//time_start = clock();
	struct  timeval start;
    struct  timeval end;
    double diff;
    gettimeofday(&start, NULL);

	// the nexus tree file (not divided, the whole)
	char * predir = getcwd(NULL, 0);
    int n = strlen(predir) + strlen(NEXUSTREE) + 2;
    char filename[n];
    i = 0;
    while((filename[i]=predir[i++])!='\0');
    strcat(filename, "/");
    strcat(filename, NEXUSTREE);
	FILE * fp;
	if((fp = fopen(filename, "r")) == NULL)
	{
		printf("Errors when opening the nexus trees!\n");
		return(-1);
	}
	fclose(fp);  // I just use the above several lines to detect whether the nexus tree file exists

	// the filehandle to report the IBD results
	FILE * fp_result;
    n = strlen(predir) + 8 + strlen(NEXUSTREE) + 1;
    char fileresult[n];
    fileresult[0] = '\0';
    strcat(fileresult, predir);
    free(predir);
    strcat(fileresult, "/result_");
    strcat(fileresult, NEXUSTREE);
    if((fp_result = fopen(fileresult, "w")) == NULL){
        printf("Cannot open the file for saving IBD results!\n");
        exit(1);
    }

    // get the TREELEN and the total tree number, and if we want several threads, we will divide the whole file
    // into several small files; and get the number of chromosomes - SAMPLE
    tree_preprocessing(filename);

    // multi-threads begin from here
    if(THREADS == 1)
    {
    	// initialize the locker
    	pthread_mutex_init(&mut, NULL);
	    int temp;
	    memset(&thread, 0, sizeof(thread));

	    // start the two threads
	    if((temp = pthread_create(&thread[0], NULL, (void *)IBD_report, (void *)fp_result)) != 0)
			printf("'IBD_report' thread not successful!\n");
	    else
			printf("'IBD_report' thread successful!\n");

		package p;
	    getpackage(&p, 1);
	    if((temp = pthread_create(&thread[1], NULL, (void *)tree_processing, (void *)&p)) != 0)
			printf("'Tree_processing' thread not successful!\n");
	    else
			printf("'Tree_processing' thread successful!\n");

		// waiting for the threads to terminate, then combine if needed
	    if(thread[1] != 0)
	    {
	        pthread_join(thread[1], NULL);
	        printf("'Tree_processing' thread over!\n");
	    }
	    else
	    	printf("'Tree_processing' thread over!\n");

	    freepackage(&p);

	    if(thread[0] != 0)
	    {
	    	IBD_report_over = 1;
	        pthread_join(thread[0], NULL);
	        printf("'IBD_report' thread over!\n");
	    }
	    else
	    	printf("'IBD_report' thread over!\n");
    }
    else  // there are several threads and boundary, we should check the boundary and combine if possible
    {
    	// initialize the locker
    	pthread_mutex_init(&mut, NULL);
	    int temp;
	    memset(&thread, 0, sizeof(thread));

	    // start the two threads
	    if((temp = pthread_create(&thread[0], NULL, (void *)IBD_report, (void *)fp_result)) != 0)
			printf("'IBD_report' thread not successful!\n");
	    else
			printf("'IBD_report' thread successful!\n");

		package p[THREADS];
		for(i=0; i<THREADS; i++)
		{
			FINISH_TABLE[i] = 0;
		    getpackage(&p[i], i+1);
		    if((temp = pthread_create(&thread[i+1], NULL, (void *)tree_processing, (void *)&p[i])) != 0)
				printf("'Tree_processing' thread #%d not successful!\n", i+1);
		    else
				printf("'Tree_processing' thread #%d successful!\n", i+1);
		}

		// waiting for the threads to terminate, then combine if needed
		// we can't use this method to combine, because there are some segments longer than the range of one thread
		/*
		int count = 0;  // we need get to (THREADS - 1) for count
		int i = 0;
		while(1)
		{
			if((FINISH_TABLE[i] == 1) && (FINISH_TABLE[i+1] != 0))  // have some problems here!!!
			{
				combine(&p[i], &p[i+1]);
				printf("Finish combining the boundary for thread #%d and thread #%d.\n", i+1, i+2);
				FINISH_TABLE[i] = 2;
				count++;
			}
			if(count == THREADS - 1)
				break;
			i++;
			if(i > (THREADS - 2))
				i = 0;
		}
		*/
		// to monitor whether or not we have finished all the trees
		int count = 0;
		i = 0;
		while(1)
		{
			if(FINISH_TABLE[i] == 1)
			{
				FINISH_TABLE[i] = 2;
				count++;
			}
			if(count == THREADS)
				break;
			i++;
			if(i > (THREADS - 1))
				i = 0;
		}

		// process the boundary, say, merge all the parts
		combine_all(p);

	    for(i=0; i<THREADS; i++)
	    {
	    	freepackage(&p[i]);
	    }

	    //========================== to remove all the temporary tree files ===========================
	    char * predir = getcwd(NULL, 0);
	    int n = strlen(predir) + strlen(NEXUSTREE) + 2 + 6;  // + '_tempx'
	    char filename[n];
	    i = 0;
	    while((filename[i]=predir[i++])!='\0');
	    strcat(filename, "/");
	    strcat(filename, NEXUSTREE);
	    strcat(filename, "_temp");
	    int filename_len = strlen(filename);
	    for(i=0; i<THREADS; i++)
	    {
	    	filename[filename_len] = (char)(i+49);
	    	filename[filename_len + 1] = '\0';
			remove(filename);  // success -> 0; fail -> -1; reason in errno
			printf("Temporary tree file '%s' has been removed.\n", filename);
	    }
	    //=============================================================================================

	    if(thread[0] != 0)
	    {
	    	IBD_report_over = 1;
	        pthread_join(thread[0], NULL);
	        printf("'IBD_report' thread over!\n");
	    }
	    else
	    	printf("'IBD_report' thread over!\n");
    }

	//time_finish = clock();
	gettimeofday(&end, NULL);
	fclose(fp_result);

	printf("Work done for tree %s!...\n", filename);
	//printf("Time used to get IBD is %f seconds.\n", (double)(time_finish - time_start)/CLOCKS_PER_SEC);
	diff = (double)(end.tv_sec-start.tv_sec)+ (double)(end.tv_usec-start.tv_usec)/1000000;
    //printf("Time used to split the original tree file is %f seconds.\n", (double)(time_finish1 - time_start)/CLOCKS_PER_SEC);
    printf("Time used totally is %f seconds.\n", diff);
	return(0);
}