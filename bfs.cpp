#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cilk.h>
#include <cilk_mutex.h>
#include <iostream>
#include "cilkview.h"
#include <sys/time.h>

cilk::cilkview cv;
using namespace std;
double elapsed_seconds()
{
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
}



int r=24;
int co = 50;
double time1, time2, elapsed;

extern "C++" {
	class node	{
			public:
				node() : v(-1), right(NULL), left(NULL) {};
				node(int vertex) : v(vertex), right(NULL), left(NULL) {};
			int v;
			node* right;
			node* left;
		
		};


	class Bag
		{
			
					
			public:
				node** list;
				int num;
				Bag()
				{
					list = new node*[r]();
					num=0;
					}
				
				
				
				node* pennant_union(node* x, node* y)
				{
					y->right = x->left;
					
					x->left = y;
					return x;
					
				}
				node* pennant_split(node* x)
				{
					node* y = x->left;
					x->left = y->right;
					y->right = NULL;
					return y;
					
				}
				
				void full_adder(node* &s1, node* &s2, node* &y)
				{
					int first;
					int second;
					int third;
					first = y!=NULL;
					second = s2!=NULL;
					third = s1!=NULL;
					
					int result =  2*2*third + 2*second + first;
					switch (result)
					{
						case 0:
							break;
							
						case 1:
							s1 = y;
							y = NULL;
							break;
						case 2:
							s1 = s2;
							break;
							
						case 3:
							pennant_union(y, s2);
							break;
							
						case 4:
							break;
							
						case 5:
							pennant_union(y, s1);
							s1 = NULL;
							break;
						
						case 6:
							pennant_union(s1,s2);
							y = s1;
							s1=NULL;
							break;
						
						case 7:
							pennant_union(y,s2);
							break;
						
					}
					
				}
				
				void Bag_insert(node* x)
				{
											

					int k=0;
					
					while(list[k]!=NULL)
					{
						x = pennant_union(x,list[k]);
						
						list[k++] = NULL ;
					}
					list[k] = x;
					
					num++;
					
				};
				
				void Bag_union(Bag s2)
				{
					node* y = NULL;
					for(int k=0; k<r; k++)
					{
						
						full_adder(list[k],s2.list[k],y);
						
					}
					num+=s2.num;
				}
				
				Bag Bag_split()
				{
					Bag s2 = Bag();
					node* y = list[0];
					list[0] = NULL;
					for(int k=0; k<r; k++)
					{
						if(list[k]!=NULL)
						{
							s2.list[k-1] = pennant_split(list[k]);
							list[k-1] = list[k];
							list[k] = NULL;
							
						}
						
					}
					
					if(y!=NULL)
						{Bag_insert(y);num--;}
						
					s2.num = num/2;
					num = num - s2.num;
					return s2;
					
				}
				
				void printBag()
				{
					for(int k=0; k<r; k++)
					{
						printTree(list[k]);
						
					}
					
				}
				
				void printTree(node* root)
				{
					if(root==NULL)
						return;
					if(root->left!=NULL)
					printTree(root->left);
					
					cout<<root->v<<endl;
					
					if(root->right!=NULL)
					printTree(root->right);	
				}
				
				node* pop()
				{
					int k=0;
					while(list[k]==NULL&&k<r)
						k++;
					if(k==r)
						return NULL;
					node* popNode =pop_back(list[k]);
					if (popNode == list[k])
						list[k] = NULL;
					return popNode;
					
				}
				
				node* pop_back(node* root)
				{
					if(root==NULL)
						return NULL;
						
					if(root->left==NULL && root->right ==NULL)	{
						return root;
						
					}
					node* popNode;
					popNode = pop_back(root->left);
					if(popNode!=NULL )
					{
						if(root->left == popNode)
						root->left=NULL;
						return popNode;
					}
					else
						popNode = pop_back(root->right);
						if(popNode!=NULL )
						{
						if(root->right == popNode)
						root->right=NULL;
						return popNode;
						}
						
					
					
				}

		};
		
		

	class Bag_reducer	{
			public:
			struct Monoid: cilk::monoid_base<Bag>	{
			
			static void reduce (Bag* left, Bag* right) {
				
				left->Bag_union(*right);
				
				}

			};
			private:
			cilk::reducer<Monoid> imp_;
			public:
			Bag_reducer() : imp_()	{}
			Bag_reducer(Bag &b) : imp_(b) {}
			//Bag_reducer(const Bag_reducer& other) : imp_(other.getValue()) {}

			Bag getValue () const
			{
				return imp_.view();
				
			}

			Bag get_value()
			{
				return imp_.view();
				
			}
			
			void Bag_insert(node* x)
			{
				Bag &b = imp_.view();
				b.Bag_insert(x);
				
			}
			
			void Bag_union(Bag s2)
			{
				Bag &s1 = imp_.view();
				s1.Bag_union(s2);
				
			}
			
			Bag Bag_split()
			{
				Bag &s1 = imp_.view();
				Bag s2 = s1.Bag_split();
				

				
				return s2;
			}
			
			int get_size()
			{
				Bag &s1 = imp_.view();
				return s1.num;
				
			}
			
			node* pop()
			{
				Bag &s1 = imp_.view();
				return s1.pop();
				
			}
			void printBag()
			{
				Bag &s1 = imp_.view();
				s1.printBag();
				
			}
			

			};
				
	
	


}

	

	typedef struct graphstruct { // A graph in compressed-adjacency-list (CSR) form
		  int nv;            // number of vertices
		  int ne;            // number of edges
		  int *nbr;          // array of neighbors of all vertices
		  int *firstnbr;     // index in nbr[] of first neighbor of each vtx
		  //node *nodeList;
	} graph;


	int read_edge_list (int **tailp, int **headp) {
		  unsigned long int max_edges = 100000000;
		  int nedges, nr, t, h;
		  *tailp = (int *) calloc(max_edges, sizeof(int));
		  *headp = (int *) calloc(max_edges, sizeof(int));
		  nedges = 0;
		  nr = scanf("%i %i",&t,&h);
		  while (nr == 2) {
			if (nedges >= max_edges) {
			  printf("Limit of %d edges exceeded.\n",max_edges);
			  exit(1);
			}
			(*tailp)[nedges] = t;
			(*headp)[nedges++] = h;
			nr = scanf("%i %i",&t,&h);
		  }
		  return nedges;
		}


	graph * graph_from_edge_list (int *tail, int* head, int nedges) {
		  graph *G;
		  int i, e, v, maxv;
		  G = (graph *) calloc(1, sizeof(graph));
		  G->ne = nedges;
		  maxv = 0;

		  // count vertices
		  for (e = 0; e < G->ne; e++) {
			if (tail[e] > maxv) maxv = tail[e];
			if (head[e] > maxv) maxv = head[e];
		  }
		  G->nv = maxv+1;
		  G->nbr = (int *) calloc(G->ne, sizeof(int));
		  G->firstnbr = (int *) calloc(G->nv+1, sizeof(int));

		  // count neighbors of vertex v in firstnbr[v+1],
		  for (e = 0; e < G->ne; e++) G->firstnbr[tail[e]+1]++;

		  // cumulative sum of neighbors gives firstnbr[] values
		  for (v = 0; v < G->nv; v++) G->firstnbr[v+1] += G->firstnbr[v];

		  for (e = 0; e < G->ne; e++) {
			i = G->firstnbr[tail[e]]++;
			G->nbr[i] = head[e];
		  }
		  for (v = G->nv; v > 0; v--) G->firstnbr[v] = G->firstnbr[v-1];
		  G->firstnbr[0] = 0;
		  return G;
		}


	void print_CSR_graph (graph *G) {
		  int vlimit = 20;
		  int elimit = 50;
		  int e,v;
		  printf("\nGraph has %d vertices and %d edges.\n",G->nv,G->ne);
		  printf("firstnbr =");
		  if (G->nv < vlimit) vlimit = G->nv;
		  for (v = 0; v <= vlimit; v++) printf(" %d",G->firstnbr[v]);
		  if (G->nv > vlimit) printf(" ...");
		  printf("\n");
		  printf("nbr =");
		  if (G->ne < elimit) elimit = G->ne;
		  for (e = 0; e < elimit; e++) printf(" %d",G->nbr[e]);
		  if (G->ne > elimit) printf(" ...");
		  printf("\n\n");
		}



	void process_Layer(Bag s1, Bag_reducer &s2, int d, graph *G, int* &level, int* &levelsize, int* &parent,  cilk::mutex &m)
	{

		if(s1.num<co)	{

			for(int i=0; i< s1.num; i++)
			{
				node* vertex = s1.pop();

				int v = vertex->v;
				
				cilk_for(int e = G->firstnbr[v]; e < G->firstnbr[v+1]; e++)	{	

					int w = G->nbr[e]; 
					if (level[w] == -1) {
					bool set=false;

						if (level[w] == -1) {
							level[w] = d+1;
							parent[w] = v;
							set = true;

							
					}

					if(set)	{
					node* newNode = new node(w);
					s2.Bag_insert(newNode);
						}
					
					}
					
				}
					
					
				}

			
		
		return;
		
		}
		

		Bag newBag = s1.Bag_split();

		cilk_spawn process_Layer(newBag, s2, d, G, level, levelsize, parent, m);
		process_Layer(s1, s2, d, G, level, levelsize, parent, m);
		cilk_sync;
		
	
	
}


	void bfs (int s, graph *G, int **levelp, int *nlevelsp, 
				 int **levelsizep, int **parentp) {
					 
		  int *level, *levelsize, *parent;
		  int thislevel;
		  int *queue, back, front;
		  int i, v, w, e;
		  level = *levelp = (int *) calloc(G->nv, sizeof(int));
		  levelsize = *levelsizep = (int *) calloc(G->nv, sizeof(int));
		  parent = *parentp = (int *) calloc(G->nv, sizeof(int));

		
		  int blks =(G->nv-1)/co +1;
		  cilk_for(int p=0; p<blks; p++)	{
			    int len =0;
				int acc = (p+1)*co;
				if(acc>G->nv)
				len=G->nv-acc+co;
				else
				len=co;
		  for (int j = 0; j < len; j++) 	{
				level[p*co+j]  = -1;
				parent[p*co+j] = -1;

				}
				
				
			}
			

		cilk::mutex m;
		  level[s] = 0;
		  levelsize[0] = 1;
		  Bag inBag;
		  node* startPoint = new node(s);
		  inBag.Bag_insert(startPoint);
		  int d=0;
		  time1 = elapsed_seconds();
		  while(inBag.num!=0){
			  Bag_reducer outBag;
			  
			  process_Layer(inBag, outBag, d, G, level, levelsize, parent, m);
			  
			  d = d+1;
			  inBag = outBag.get_value();
			

		  }
		  time2 = elapsed_seconds();
		 
		  *nlevelsp = d;
		}


	int cilk_main (int argc, char* argv[]) {
		int p =cilk::current_worker_count();
		#pragma cilk grainsize = 1;

		  graph *G;
		  int *level, *levelsize, *parent;
		  int *levelsize2;
		  int *tail, *head;
		  int nedges;
		  int nlevels;
		  int startvtx;
		  int i, v, reached;
		startvtx = atoi (argv[1]);
		 /* if (argc == 2) {
			startvtx = atoi (argv[1]);
		  } else {
			printf("usage:   bfstest <startvtx> < <edgelistfile>\n");
			printf("example: cat sample.txt | ./bfstest 1\n");
			exit(1);
		  }*/
		  nedges = read_edge_list (&tail, &head);
		  G = graph_from_edge_list (tail, head, nedges);
		  free(tail);
		  free(head);
		  print_CSR_graph (G);
		  int total = G->nv;
		  int step=0;
		  while(total>0)
		  {
			  total = total>>1;
			  ++step;
			  
			  
		  }
		  r = step+1;
		  printf("Starting vertex for BFS is %d.\n\n",startvtx);
		   //co = 50< G->nv/(p)? 50: G->nv/(p);

			//if(co==0)
			//co=1;
		  co = 128;
		  cv.start();
			
		  bfs (startvtx, G, &level, &nlevels, &levelsize, &parent);

		  cv.stop();
		  elapsed = time2 - time1;
		  cv.dump("PBFS");
		  levelsize2 = (int *) calloc(G->nv, sizeof(int));
		  for( i =0; i< G->nv; i++)	{
			  if(level[i]!=-1)
				levelsize2[level[i]]++;
		
			}
		  reached = 0;
		  for (i = 0; i < nlevels; i++) reached += levelsize2[i];
		  printf("Breadth-first search from vertex %d reached %d levels and %d vertices.\n",
			startvtx, nlevels, reached);
		  for (i = 0; i < nlevels; i++) printf("level %d vertices: %d\n", i, levelsize2[i]);
		  if (G->nv < 20) {
			printf("\n  vertex parent  level\n");
			for (v = 0; v < G->nv; v++) printf("%6d%7d%7d\n", v, parent[v], level[v]);
		  }
		  printf("\n");
		  
		  
	 printf("current worker count %d\n",p);
	 printf("runtime : %e \n",elapsed);
	 printf("coarse: %d \n", co);
		 
		}




