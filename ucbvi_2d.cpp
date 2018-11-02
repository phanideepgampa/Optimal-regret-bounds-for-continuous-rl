
#include <math.h>
#include <bits/stdc++.h>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <cstdlib>
#include <random>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include "serialize_tuple.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace std;
using namespace boost::archive;
const long double PI=3.141592653589793238462643383279502884197169399375105820974944592307816406286;
int n=5,A=2,K=1200000,d=2;
int sub=n/5;
double H=20;
double T= K*H;
double S =pow(n,d),delta=0.100,L =log((5*(S)*A*T)/delta), b_temp= ((pow(10,4))*(pow(H,3))*(pow(S,2))*A*(pow(L,2)));
const string path="regret_2d/5_5_20";
double x=0,y=0;
int i1,i2=0;
mt19937 gen(1234),gen2(123);
uniform_real_distribution<> dis(-(PI/(2*sqrt(2)))+0.00000001, ((19.0*PI)/(2*sqrt(2)))-0.000000001);
uniform_real_distribution<> dis2((PI/(2*sqrt(2)))+0.00000001,((21.0*PI)/(2*sqrt(2)))-0.000000001);//root(2)*pi is the spacing
vector< vector<double> > r((n*n*n), vector<double>(2) );//(x,v)
int ids[300];
unordered_map < pair< int , int > , unordered_set< int > ,boost:: hash<pair<int,int> > > prob_helper; //use id  returns vector of y ids
unordered_map < pair< int , int >, double , boost :: hash<pair<int,int> > >  kappa;// use id
unordered_map < pair<int , int> , double ,boost:: hash<pair<int,int> > > n_prime ;  //use id
unordered_map < tuple<int,int, int> ,double ,boost:: hash<tuple<int,int, int>  > > probcount;//use id
unordered_map <  tuple<int,int ,int>, double ,  boost ::hash<tuple<int,int, int> >  > probability; //use id
unordered_map < pair<int,int>, double, boost:: hash<pair<int,int> > > V;
unordered_map < pair<int,int>, bool, boost:: hash<pair<int,int> > > check_v;
float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r1 = random * diff;
    return a + r1;
}

double Q[1000][2][102]; 
void take_action()
{
	x=dis(gen);
	y=dis2(gen2);
}
int intervals_to_id(int i,int j){
	int p=(i+j)*(i+j+1);
	p>>=1;
	
	 p+=j;
	                            //cantor pairing function;
	return p;
}

void state_to_interval()
{  
    long double x_origin=-(PI/(2*sqrt(2))),y_origin=(PI/(2*sqrt(2))),diff=(PI*(sqrt(2)));
	i1=ceil(((x-x_origin)*sub)/diff);
	i2=ceil(((y-y_origin)*sub)/diff);
	
	
}
long double fn1(double p,double q)
{
	 return (sin((p-q)/sqrt(2))+cos((p+q)/sqrt(2)));
}
long double fn2(double p,double q)
{
	 return (-sin((p-q)/sqrt(2))-cos((p+q)/sqrt(2)));
}
void init_reward(){
			
	double max1=-10000,min1=10000000,max2=-10000,min2=100000000;
  int m=1;
   double x_origin=-(PI/(2*sqrt(2))),y_origin=(PI/(2*sqrt(2))),diff=(PI*(sqrt(2)));
   	double left_x=x_origin,left_y=y_origin,add=diff/double(sub),inte1=0,inte2=0;
	for (int i=1;i<=n;i++)
	{  left_y=y_origin;
		for(int j=1;j<=n;j++)
		 {
		 	
	        double left_sx=left_x,left_sy=left_y,add_s=add/double(n);
	   		 inte1=0,inte2=0;
	   
	     for(int dx=1;dx<=n;dx++)
	     
	     {  left_sy=left_y;
	     	
	     	 for(int dy=1;dy<=n;dy++)
	     	 {
	     	 	
	     	 	if(fn1(left_sx,left_sy)<fn1(left_sx+add_s,left_sy+add_s))
	     	 	 {
	     	 	 	inte1+=(fn1(left_sx,left_sy));
				   }
				else{
					inte1+=fn1(left_sx+add_s,left_sy+add_s);
					
				}
	     	 	
	     	 	if(fn2(left_sx,left_sy)<fn2(left_sx+add_s,left_sy+add_s))
	     	 	 {
	     	 	 	inte2+=(fn2(left_sx,left_sy));
				   }
				else{
					inte2+=fn2(left_sx+add_s,left_sy+add_s);
					
				}
				
	     	 	left_sy+=add_s;
			  }
	     	 left_sx+=add_s;
	     	
		 }
	       
	       
	       
	     inte1/=double(n*n);
		 inte2/=double(n*n);
		 max1=max(max1,inte1);
		 max2=max(max2,inte2);
		 min1=min(min1,inte1);
		 min2=min(min2,inte2);
		 int id=intervals_to_id(i,j);
		 r[id][0]=inte1;
		 r[id][1]=inte2;
	      //cout<<left_x<<" "<<left_y<<" "<<r[id][0]<<" "<<r[id][1]<<endl; 
		 	left_y+=add;
		 }
		 
		left_x+=add;
		
		 
		}
		
		for(int i=1;i<=n;i++)
		{
			for(int j=1;j<=n;j++)
		   {  int id=intervals_to_id(i,j);
			 r[id][0]=((r[id][0]-min1)/(max1-min1));
			 r[id][1]=((r[id][1]-min2)/(max2-min2));
			cout<<r[id][0]<<" "<<r[id][1]<<endl;
		   }
		}
	
	
	int n_3=n*n*n;
	for(int i=0;i<n_3;i++)
	 for(int j=0;j<2;j++)
	  for(int k=0;k<102;k++)
	    Q[i][j][k]=H;
	
	
}
/*BONUS*/

pair <double,double> bonus_BF(int x_id,int a,int h)
{   pair <int,int> x_a= make_pair(x_id,a);
	
	 double mean=0,rightterm=0;
	 for(auto it =prob_helper[x_a].begin();it!= prob_helper[x_a].end();++it )
	 {  int y_id=*it;
	 	 if(check_v[make_pair(y_id,h+1)]){
	 	 	   
	 	 	  mean+=(V[make_pair(y_id,h+1)]*probability[make_tuple(y_id,x_id,a)]);
	 	 	
	 	 	
		  }
		  else if(h!=H){
		  	
		  	mean+=(double(H)*probability[make_tuple(y_id,x_id,a)]);
		  	
		  	
		  }
	 	
	 	
	 	if(n_prime[make_pair(y_id,h+1)])
	 	  {
	 	  	
	 	  	 rightterm+=(probability[make_tuple(y_id,x_id,a)]*min(pow(H,2),(b_temp/n_prime[make_pair(y_id,h+1)])));
	 	  	
		   }
	 	else{
	 		rightterm+=(probability[make_tuple(y_id,x_id,a)]*pow(H,2));
	 		
		 }
	 	
	 	
	 	
	 }
	 
	 double variance=0;
	 
	 for( auto it =prob_helper[x_a].begin();it!= prob_helper[x_a].end();++it )
	 {  
	   int y_id=*it;
	 	 if(check_v[make_pair(y_id,h+1)]){
	 	 	   
	 	 	  variance+=((V[make_pair(y_id,h+1)]-mean)*(V[make_pair(y_id,h+1)]-mean)*probability[make_tuple(y_id,x_id,a)]);
	 	 	
	 	 	
		  }
		  else if(h!=H){
		  	
		  	variance+=((H-mean)*(H-mean)*probability[make_tuple(y_id,x_id,a)]);
		  	
		  	
		  }
		  else{
		  	
		  	variance+=((mean)*(mean)*probability[make_tuple(y_id,x_id,a)]);
		  	
		  }
	 	
	 	
	 	
	 }
	 
	 
	 double bonus=sqrt((8*L*variance)/(kappa[x_a]))+((14*H*L)/(3*kappa[x_a]))+sqrt((8*rightterm)/(kappa[x_a]));
	 
	 return make_pair(bonus,mean);
	 
	
	
	
}

void output_Q(){
	
	ofstream output;
	
	output.open(path+"/Q.txt");
	int n_3=n*n*n;
	for(int i=0;i<n_3;i++)
	 {
	 for(int j=0;j<2;j++){
	 
	  for(int k=0;k<102;k++)
	   {
		      output<<Q[i][j][k]<<":";
	}
	output<<"?";
	
}
output<<endl;

}
	output.close();
	
}


void UCBVI()
{  ofstream file1,file2;
   init_reward();
	
	file1.open(path+"/file1.txt");
	file2.open(path+"/file2.txt");
	cout<<"reward done";
	 for(int episode=0;episode<K;episode++)
	 {  
	    auto start = std::chrono::high_resolution_clock::now();

	    
		    double av_bonus=0,av_mean=0;

		 
		 for(int horizon=H;horizon>0;horizon--)
		 {
		 	for ( auto it = kappa.begin(); it != kappa.end(); ++it )
		 	{    int x_id=it->first.first;
		 	   int a=it->first.second;
		 		  pair<double,double> bonus= bonus_BF(x_id,a,horizon);
		 		  //cout<<"---"<<r_id[x_id]+bonus.first+bonus.second<<endl;
		 		av_bonus+=bonus.first;
		 		av_mean+=bonus.second;
		 		Q[x_id][a][horizon]=min(Q[x_id][a][horizon],min(double(H),r[x_id][a]+bonus.first+bonus.second));
		 		//cout<<" ---"<<Q[x_id][a][horizon]<<endl;
		 		if((Q[x_id][a][horizon]!=H)&&(Q[x_id][1-a][horizon]!=H))
				 
		 		 { 
				    V[make_pair(x_id,horizon)]=max(Q[x_id][a][horizon],Q[x_id][1-a][horizon]);
				    check_v[make_pair(x_id,horizon)]=1;
				}
			 }
		 	
		 	}
		av_bonus/=(H*kappa.size());
		
		av_mean/=(H*kappa.size());
			 
			int count=0,a_left=0,a_right=1;
			double avreward=0,av_x=0,acc=0;
		
			take_action();
			double res1=x,res2=y,v_opt=0,v_pol=0,v_opt2=0,v_pol2=0;
			int his=(rand()%10);
		for(int exp=0;exp<10;exp++)	
	 	{
		    x=res1;
		    y=res2;
		 for(int horizon=1;horizon<=H;horizon++)
	 	{
	 		state_to_interval();
	 		int x_id=intervals_to_id(i1,i2); 
	 		int action;
	 		double decider=RandomFloat(0,1);
	 		 
	 		   if((Q[x_id][0][horizon]==H)&&(Q[x_id][1][horizon]==H))
	 		    {  
				   action=decider>0.5?0:1;
				   
			}
	 		else
          		action=Q[x_id][0][horizon]>Q[x_id][1][horizon]? 0:1;
	 		
			 
			 
			if(action)
				v_pol+=((2+fn2(x,y))/4.0);
			else
			  v_pol+=((2+fn1(x,y))/4.0);
			
	 	   
	 	    if(fn1(x,y)>fn2(x,y))
	 	      v_opt+=((2+fn1(x,y))/4.0);
	 	    else
	 	      v_opt+=((2+fn2(x,y))/4.0);
	 		
	 		
	 		 
			  take_action();
			 state_to_interval();
			 int y_id;  
			  
	 		 y_id=intervals_to_id(i1,i2);
	 		      	
		
			
	 	   //file1<<x_id<<" "<<y_id<<" "<<action<<endl;
		 if(exp==his)
	 	   {
			 
	 		kappa[make_pair(x_id,action)]++;
	 		n_prime[make_pair(y_id,horizon)]++;
	 		probcount[make_tuple(y_id,x_id,action)]++;
	 		
	 		probability[make_tuple(y_id,x_id,action)]=(probcount[make_tuple(y_id,x_id,action)])/(kappa[make_pair(x_id,action)]);
	 		prob_helper[make_pair(x_id,action)].insert(y_id);
	 		
	      }
	 		
	 	 
	 	   
		 }
	 	

	 	
	 }
	 v_opt/=(10.0);
	 v_pol/=(10.0);
	 //v_opt2/=(10.0);
	 //v_pol2/=(10.0);
	 double regret=v_opt-v_pol;
	 
	 	auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
	  file1<<"------- "<<episode<<"  "<<regret<<endl;
	 	file2<<"--- "<<episode<<" "<<V.size()<<" "<< kappa.size()<<" "<<"bonus "<<av_bonus<<" mean "<<av_mean<<" "<<elapsed.count()<<"sec"<<endl;
	 	V.clear();
	 	check_v.clear();
	 	
	 	//output_Q();
	 	
	
}
	file1.close();
	file2.close();
	
}



void save_history(){
	ofstream h1(path+"/prob_helper"),h2(path+"/kappa"),h3(path+"/n_prime"),h4(path+"/probcount"),h5(path+"/probability");
	
	text_oarchive o1(h1),o2(h2),o3(h3),o4(h4),o5(h5);
	o1<<prob_helper;
	o2<<kappa;
	o3<<n_prime;
	o4<<probcount;
	o5<<probability;
	h1.close();
	h2.close();
	h3.close();
	h4.close();
	h5.close();
}
void load_history(){
	ifstream l1(path+"/prob_helper"),l2(path+"/kappa"),l3(path+"/n_prime"),l4(path+"/probcount"),l5(path+"/probability");
	
	text_iarchive i1(l1),i2(l2),i3(l3),i4(l4),i5(l5);
	i1>>prob_helper;
	i2>>kappa;
	i3>>n_prime;
	i4>>probcount;
	i5>>probability;
	l1.close();
	l2.close();
	l3.close();
	l4.close();
	l5.close();

	
}




int main()
{
	srand (static_cast <unsigned> (1234));
	//init_reward();
  

// load_history();
UCBVI();
 
output_Q();
save_history();
 //load_history();
 //init_reward();
 
  //for(auto it=probability.begin();it!=probability.end();++it)
    //{
    //	cout<<it->second<<endl;
	//}

	
}




