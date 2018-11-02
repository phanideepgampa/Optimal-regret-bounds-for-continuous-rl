#include <math.h>


#define GRAVITY 9.8
#define MASSCART 1.0
#define LENGTH 0.5		  
//#define FORCE_MAG 1.5
#define TAU 0.02		  /* seconds between state updates */
#define x_threshold 3.00000
//#define v_threshold 1.50000 


#include <bits/stdc++.h>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <cstdlib>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include "serialize_tuple.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace std;
using namespace boost::archive;

int n=6,A=2,K=1000000,d=2;
double H=100;
double T= K*H,init=0.04;
double S =pow(n,d),delta=0.100,L =log((5*(S)*A*T)/delta), b_temp= ((pow(10,4))*(pow(H,3))*(pow(S,2))*A*(pow(L,2)));
double v_threshold=(x_threshold-init)/(H*TAU);
double FORCE_MAG=(v_threshold-init)/(H*TAU);
double x=0,v=0;
int i1,i2;

vector< vector<double> > r(n+1, vector<double>(n+1) );//(x,v)
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
    float r = random * diff;
    return a + r;
}

//H+1
int intervals_to_id(int i,int j){
	int p=(i+j)*(i+j+1);
	p>>=1;
	
	 p+=j;
	                            //cantor pairing function;
	return p;
}
int solution(){
	
	
	if(-x>((v*abs(v))/(2*FORCE_MAG))) return 1;
	//	if(-x>((v*v)/(2*FORCE_MAG))) return 1;

	else return 0;
	  
	
}
int qsize=intervals_to_id(n+1,n+1);
vector< double> r_id(qsize);
double Q[300][2][102]={double(H)}; //change dim when n changes,H
int take_action( int a )
{
	
	if(x>=x_threshold|| x<=(-x_threshold)|| (v>=v_threshold)||(v<=(-v_threshold)) )
	   return -1;
	   
	double force = (a>0)? FORCE_MAG : -FORCE_MAG;
	double acceleration=(force/MASSCART);
	x+=(TAU*v);
	v+=(TAU*acceleration);
	
	
	return 1;
}
void reset(){
	//x = -x_threshold + (static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2*x_threshold))));
	//v = -v_threshold + (static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2*v_threshold))));
    x=RandomFloat(-init,init);
    v=RandomFloat(-init,init);
}


void state_to_interval()
{
	 double tempx=double((x+x_threshold)/(2*x_threshold));
	 double tempv=double((v+v_threshold)/(2*v_threshold));
	  
	 i1=ceil(tempx*n);
	 i2=ceil(tempv*n);
	
}

void init_reward(){
	double right_i,left_i,right_j,left_j,addx,addv,reward;
	double maxxv=x_threshold*v_threshold;
	double mini=10000000,maxi=-1000000;
	for (int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
			 right_i=(((i*2*x_threshold)/n)-x_threshold);
			 left_i=double((((i-1)*2*x_threshold)/n)-x_threshold);
			 right_j=double(((j*2*v_threshold)/n)-v_threshold);
			 left_j=double((((j-1)*2*v_threshold)/n)-v_threshold);
			 addx=(right_i-left_i)/n;
			 addv=(right_j-left_j)/n;
			 reward=0;
			 //cout<<left_i<<" "<<right_i<<endl;
			for(int k=0;k<n;k++)
			{ left_i+=addx;
			  left_j+=addv;
			  reward+=(maxxv-abs(left_i*left_j));
			}
				reward/=(n*n);
			//cout<<reward<<endl;
			maxi=max(maxi,reward);
			mini=min(mini,reward);
			r[i][j]=reward;
			
			
			
	}
			
			
			
	}
	
	
	for (int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
		  r[i][j]=(r[i][j]-mini)/(maxi-mini); 
		  //cout<<r[i][j]<<endl;
		  r_id[intervals_to_id(i,j)]=r[i][j]; 
		}
		
	}
	
	for(int i=0;i<300;i++)
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
	
	output.open("n_6_r/Q_1_100.txt");
	
	for(int i=0;i<300;i++)
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
	
	file1.open("n_6_r/file1_100.txt");
	file2.open("n_6_r/file2_100.txt");
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
		 		Q[x_id][a][horizon]=min(Q[x_id][a][horizon],min(double(H),r_id[x_id]+bonus.first+bonus.second));
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
			
			reset();
	 	for(int horizon=1;horizon<=H;horizon++)
	 	{
	 		state_to_interval();
	 		int x_id=intervals_to_id(i1,i2); 
	 	
	 	    
	 	    
			 
	 		int action;
	 		double decider=RandomFloat(0,1);
	 		 
	 		   if(Q[x_id][0][horizon]==H&&Q[x_id][1][horizon]==H)
	 		     action=decider>0.5?0:1;
	 		else
          		action=Q[x_id][0][horizon]>Q[x_id][1][horizon]? 0:1;
	 		
			 av_x+=x;
	 		if(action)a_right++;
	 		else a_left++;
	 		
	 		 if(action==solution())
	 		   acc++;
			 int done= take_action(action);
			 state_to_interval();
			 int y_id;
	 	
	 	        /*  
			     if(x>=x_threshold)
			       {  double low=x_threshold-((2*x_threshold)/double(n));
			       	   x=RandomFloat(low+0.0001,x_threshold-0.000001);
				   }
				 else if(x<=(-x_threshold))
				 {
				 	double hi=-x_threshold+((2*x_threshold)/double(n));
			       	   x=RandomFloat(-x_threshold+0.00001,hi-0.00001);
				 }
			  
			  if(v>=v_threshold)
			       {  double low=v_threshold-((2*v_threshold)/double(n));
			       	   v=RandomFloat(low+0.0001,v_threshold-0.000001);
				   }
				 else if(v<=(-v_threshold))
				 {
				 	double hi=-v_threshold+((2*v_threshold)/double(n));
			       	   v=RandomFloat(-v_threshold+0.00001,hi-0.00001);
				 }
			  
			  state_to_interval();
			   
			  */
			  
	 		 y_id=intervals_to_id(i1,i2);
	 		 
			     	
			 
	 	   
	 	      
	 	   
	 	    avreward+=r_id[y_id];
	 	    
	 		kappa[make_pair(x_id,action)]++;
	 		n_prime[make_pair(y_id,horizon)]++;
	 		probcount[make_tuple(y_id,x_id,action)]++;
	 		
	 		probability[make_tuple(y_id,x_id,action)]=(probcount[make_tuple(y_id,x_id,action)])/(kappa[make_pair(x_id,action)]);
	 		prob_helper[make_pair(x_id,action)].insert(y_id);
	 		
	 	if((i1==(n/2)||i1==((n+2)/2))&&(i2==(n/2)||i2==((n+2)/2)))
	 		  count++;
	 		
	 	 
	 	   
		 }
	 	
	 	av_x/=H;
	 	 avreward/=H;
	 	 acc/=H;
	 	auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
	  file1<<"------- "<<episode<<"  "<<count<<" "<<a_left<<" "<<a_right<<" "<<V.size()<<" "<<acc<<endl;
	 	file2<<"--- "<<episode<<" "<<avreward<<" "<< kappa.size()<<" "<<"bonus "<<av_bonus<<" mean "<<av_mean<<" x "<<av_x<<" "<<elapsed.count()<<"sec"<<endl;
	 	V.clear();
	 	check_v.clear();
	 	
	 	output_Q();
	 	
	 	
	 }
	
	
	file1.close();
	file2.close();
	
}



void save_history(){
	ofstream h1("n_6_r/prob_helper_100"),h2("n_6_r/kappa_100"),h3("n_6_r/n_prime_100"),h4("n_6_r/probcount_100"),h5("n_6_r/probability_100");
	
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
	ifstream l1("n_6_r/prob_helper_100"),l2("n_6_r/kappa_100"),l3("n_6_r/n_prime_100"),l4("n_6_r/probcount_100"),l5("n_6_r/probability_100");
	
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
  /* for(int i=1;i<=6;i++)
     for(int j=1;j<=6;j++)
       cout<< i<<" "<<j <<" "<<r_id[intervals_to_id(i,j)]<<endl;
  */ 
   
  //reset();
  
  //cout<<take_action(0)<<endl;
   	//cout<<i1<<" "<<i2<<endl;
 // load_history();
 cout<<v_threshold<<" "<<FORCE_MAG<<endl;
 UCBVI();
 //init_reward();
//for (int i=0;i<100;i++)
  //cout<<RandomFloat(-60.0,60.0)<<endl;	
 //cout<<r[50][50]<<" "<<r[48][49]<<" "<<r[52][51];
 
 output_Q();
 save_history();
 //load_history();
      

	
}




