#ifndef COUPLE_CELLS_SIMULATION    // To make sure I don't declare the function more than once by including the header multiple times.
#define COUPLE_CELLS_SIMULATION

#include <iostream>
#include <fstream>
#include <cmath>
#include<vector>

using namespace std;


int defaultInt[2] = {0,1};
float defaultFloat[2] = {0,1};


/*Define some parameters */
const float R_F = 8.6174e-005; //needed for the computation of the revsersal for Ca
const float T = 10;
double dt = 0.01;            //time step in ms 




float getThreshold(vector<float>& spike_times_, int length_of_spike_times);
bool evaluate_spike_train(vector<float>& spike_train_0, vector<float>& spike_train_1, 
              int length_of_spike_train_0 , int length_of_spike_train_1 ,
              int simulation_time_, 
              float min_period_threshold, float max_period_threshold,
              float min_phase_difference_threshold, float max_phase_difference_threshold,
              float min_burst_exclusion_threshold);

double get_power(double x, int y)
{

  double prod = 1.0;
  for (int kk=0;kk<y;kk++)
    {
      prod = prod*x;
    }
  return prod;
}
double Xinf(double V0, double A, double B, double C = 1)
{

  double X = C/(1.0+exp((V0+A)/B));
  return X;
}


double Xinf_dV(double V0, double A, double B, double C = 1)
{
  double X = -C*exp((V0+A)/B)/(B*get_power((1.0+exp((V0+A)/B)),2));
  return X;
}


double tauX(double V0, double A, double B, double D, double E)
{
  double X = A - B/(1+exp((V0+D)/E));
  return X;
}

double tauhNa(double V0)
{
  double X = (0.67/(1.0+exp((V0+62.9)/-10.0)))*(1.5 + 1.0/(1.0+exp((V0+34.9)/3.6)));
  return X;
}

double taumCaS(double V0)
{
  double X = 1.4 + (7.0/((exp((V0+27.0)/10.0))+(exp((V0+70.0)/-13.0))));
  return X;
}

double tauhCaS(double V0)
{
  double X = 60.0 + (150.0/((exp((V0+55)/9.0))+(exp((V0+65.0)/-16.0))));
  return X;
}



void read_stimulus(FILE *stimfile, double g_vec[8][2001])
{ 
  /*Read in stimulus times from a file provided.*/
 double x = 1.0;

 for (int ii=0;ii<2001;ii++)
   {
     for (int jj=0;jj<8;jj++)
       {
   fscanf(stimfile,"%lf \t",&x);
   g_vec[jj][ii] = x;
       }
   }
 
 return;
}






class neuron {
   
  public:
    neuron(int cell_number_, double V_, double g_syn_ , double E_syn_, bool supress_output_, int max_con_to_change[], float change_max_con_to[], int how_many_changes);
    void evolve_gates( void );
    void evolve_V(double Svar_ ,double input_current_ );
    void evolve_Ca( void );
    void calculate_time_constants( void );
    double get_DIC(int ii_);
    void weight(int ii_);
    void calculate_DIC( void );
    void set_to_infinity( void );
    double get_V( void );
    double get_Ca( void );
    void set_V( double V_ );
    double get_membrane_current( void );
    double Ca_dV( void );
    double KCainf_dV( bool Ca_ );

  private:
    
    double V;
    double I;
    double Ca; 
    double Ca_inf;
    double I_Ca;
    float Capacity;
    double gate[2][7];
    double tau[2][8];
    double tau_Ca;
    double tau_V;
    double gate_inf[2][7];
    double E[7];
    double E_L; 
    double E_syn;
    double g[7];
    double g_L;
    double g_syn;
    int p[7];
    double gnmh[7];
    double gnmh_Ca;
    double sum_g;
    double sum_Eg;
    double V_inf;
    double weight_array[4];
    double DIC[3];
    double Xinf_par[7][2][2];

    bool incoming_connection;

  


};

// class member functions

neuron::neuron(int cell_number_, double V_ , double g_syn_ = 0, double E_syn_ = 0 , bool supress_output_ = false, int max_con_to_change[] = defaultInt , float change_max_con_to[] = defaultFloat, int how_many_changes = 0){



  if (g_syn_ == 0 and E_syn_ == 0){
    // in this case there is no incoming connection
    incoming_connection = false;
  }
  else{
    incoming_connection = true;
    g_syn = g_syn_;
    E_syn = E_syn_;
  }

  // parameters taken from the liu paper
    Xinf_par[0][0][0]  = 25.5;
    Xinf_par[0][0][1]  = -5.29; // m, Na
    Xinf_par[0][1][0]  = 48.9;
    Xinf_par[0][1][1]  = 5.18; // h, Na
    Xinf_par[1][0][0]  = 27.1;
    Xinf_par[1][0][1]  = -7.2; // m, CaT
    Xinf_par[1][1][0]  = 32.1;
    Xinf_par[1][1][1]  = 5.5; // h, CaT
    Xinf_par[2][0][0]  = 33.0;
    Xinf_par[2][0][1]  = -8.1; // m, CaS
    Xinf_par[2][1][0]  = 60.0;
    Xinf_par[2][1][1]  = 6.2; // h, CaS
    Xinf_par[3][0][0]  = 27.2;
    Xinf_par[3][0][1]  = -8.7; // m, A
    Xinf_par[3][1][0]  = 56.9;
    Xinf_par[3][1][1]  = 4.9; // h, A
    Xinf_par[4][0][0]  = 28.3;
    Xinf_par[4][0][1]  = -12.6; // m, KCa
    Xinf_par[4][1][0]  = 0;
    Xinf_par[4][1][1]  = 0; // h, KCa
    Xinf_par[5][0][0]  = 12.3;
    Xinf_par[5][0][1]  = -11.8; // m, Kd
    Xinf_par[5][1][0]  = 0;
    Xinf_par[5][1][1]  = 0; // h, Kd
    Xinf_par[6][0][0]  = 70.0;
    Xinf_par[6][0][1]  = 6.0; // m, Kd
    Xinf_par[6][1][0]  = 0;
    Xinf_par[6][1][1]  = 0; // h, Kd


  // p values
  int p_Na = 3;
  int p_A = 3;
  int p_K = 4;
  int p_H = 1;
  int p_CaT = 3;
  int p_CaS = 3;
  int p_KCa = 4;

  // reversal potentials
  double E_Na = 50.0; //mV or 30?
  double E_A = -80.0; //mV
  double E_K = -80.0; //mV
  double E_H = -20.0; //mV
  double E_Ca = 120; //mV default, but will vary

  V = V_;
  Ca = 10;
  tau_Ca = 200;
  I_Ca = 0;
  E_L = -50.0; //mV
  g_L = 0.01;
  Capacity = 1; //nF



  // get input from the file, containing different maximum input conductances
  double g_vec[8][2001];
  char stim_file [999];
  FILE *stimfile;
  sprintf (stim_file, "conductances_dataset_2000.dat");
  stimfile = fopen(stim_file,"r");
  read_stimulus(stimfile,g_vec);
  fclose(stimfile);

  for (int ii = 0; ii < 7; ++ii)
  {
    g[ii] = g_vec[ii+1][cell_number_];
    if (not supress_output_){
    cout << g[ii] << " ";
    }
    for (int jj = 0; jj < how_many_changes; ++jj)
    {
      if (max_con_to_change[jj] == ii)
      {
        g[ii] = change_max_con_to[jj];
        if (not supress_output_)
        {
        cout << " changed max g for " << ii << " to: " << g[ii] << " ";
        }
      }
    }
  } 

  E[0] = E_Na; 
  E[1] = E_Ca;
  E[2] = E_Ca;
  E[3] = E_A;  
  E[4] = E_K;  
  E[5] = E_K;  
  E[6] = E_H;  
  //
  p[0] = p_Na;
  p[1] = p_CaT;
  p[2] = p_CaS;
  p[3] = p_A;
  p[4] = p_KCa;
  p[5] = p_K;
  p[6] = p_H; 


  for (int ii = 0; ii < 7; ++ii)
  { 

    gate_inf[0][ii] = Xinf(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1]);
    if (ii == 4)
    {
      gate_inf[0][ii] = Xinf(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1],Ca/(Ca+3.0));
    }
    if (ii < 4)
    {
      gate_inf[1][ii] = Xinf(V,Xinf_par[ii][1][0],Xinf_par[ii][1][1]);
    }
    else
    {
      gate_inf[1][ii] = 1.0;
    }

    gate[0][ii] =  gate_inf[0][ii];
    gate[1][ii] =  gate_inf[1][ii];
    gnmh[ii] = 0.0;
  }

  if (not supress_output_){
    cout << endl << "Neuron " << cell_number_ <<  " has been created. " <<endl;
  }
}

void neuron::calculate_time_constants( void ){

    tau[0][0] = tauX(V,1.32,1.26,120.0,-25); // m, Na
    tau[1][0] = tauhNa(V); // h, Na
    tau[0][1] = tauX(V,21.7,21.3,68.1,-20.5); // m, CaT
    tau[1][1] = tauX(V,105.0,89.8,55.0,-16.9); // h, CaT
    tau[0][2] = taumCaS(V); // m, CaS
    tau[1][2] = tauhCaS(V); // h, CaS
    tau[0][3] = tauX(V,11.6,10.4,32.9,-15.2); // m, A
    tau[1][3] = tauX(V,38.6,29.2,38.9,-26.5); // h, A
    tau[0][4] = tauX(V,90.3,75.1,46.0,-22.7); // m, KCa
    tau[1][4] = 0.0; // h, kCa
    tau[0][5] = tauX(V,7.2,6.4,28.3,-19.2); // m, Kd
    tau[1][5] = 0.0; // h, Kd 
    tau[0][6] = tauX(V,272.0,-1499.0,42.2,-8.73); // m, H
    tau[1][6] = 0.0; // h, H
    tau[0][7] = tau_Ca;
    tau[1][7] = 0.0;



}

void neuron::evolve_gates( void ){
  for (int ii = 0; ii < 7; ++ii)
    {
      gate_inf[0][ii] = Xinf(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1]);
      if (ii == 4)
      {
        gate_inf[0][ii] = Xinf(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1],Ca/(Ca+3.0));
      }
      if (ii < 4)
      {
        gate_inf[1][ii] = Xinf(V,Xinf_par[ii][1][0],Xinf_par[ii][1][1]);
      }

    }

  for (int ii = 0; ii < 7 ; ++ii)
     {
       gate[0][ii] = gate[0][ii] + (1-exp(-dt/tau[0][ii]))*(gate_inf[0][ii] - gate[0][ii]);
       if ( ii < 4)
       {
        gate[1][ii] = gate[1][ii] + (1-exp(-dt/tau[1][ii]))*(gate_inf[1][ii] - gate[1][ii]);
       }
       
     }

}

void neuron::evolve_Ca( void ){

  // integrate Ca dynamics 
  gnmh_Ca = 0;
  for (int ii = 1; ii < 3; ++ii)
     {
       gnmh_Ca += g[ii]*get_power(gate[0][ii],p[ii])*gate[1][ii];
     }
  I_Ca = gnmh_Ca*(V-E[1]);
  Ca_inf = 0.05 - 0.94*I_Ca; 
  Ca = Ca + (1-exp(-dt/tau_Ca))*(Ca_inf - Ca);
  E[1] = E[2] = 500.0*R_F*(T + 273.15)*log(3000.0/Ca);
   
}

void neuron::evolve_V( double Svar = 0, double input_current_ = 0 ){
  if (incoming_connection)
{  sum_g = g_L + g_syn * Svar;
  sum_Eg = g_L*E_L + g_syn * Svar * E_syn ;

}
else 
{
  sum_g = g_L;
  sum_Eg = g_L*E_L;
}
   for (int ii = 0; ii < 7; ++ii)
     {
       gnmh[ii] = g[ii]*get_power(gate[0][ii],p[ii])*gate[1][ii];
       sum_g += gnmh[ii];
       sum_Eg += gnmh[ii]*E[ii];

     }
   
   tau_V  = Capacity/sum_g;  //Membrane time constant.
   // V_inf is steady-state voltage.
   V_inf = (sum_Eg - input_current_)/sum_g ;
   
   //evolve the voltage using exponential Euler
   V = V_inf + (V-V_inf)*exp(-dt/tau_V);

}

void neuron::set_to_infinity( void ){


  // this section is needed to compute the desired g_values for the calcuim current
     for (int ii = 1; ii < 3; ++ii)
    { 
        gate_inf[0][ii] = Xinf(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1]);
        gate[0][ii] = gate_inf[0][ii];
        gate_inf[1][ii] = Xinf(V,Xinf_par[ii][1][0],Xinf_par[ii][1][1]);
        gate[1][ii] = gate_inf[1][ii];
    }

  // calculating the clacium current/concentration
    gnmh_Ca = 0;
    for (int ii = 1; ii < 3; ++ii)
       {
         gnmh_Ca += g[ii]*get_power(gate[0][ii],p[ii])*gate[1][ii];
       }
    I_Ca = gnmh_Ca*(V-E[1]);
    Ca_inf = 0.05 - 0.94*I_Ca ; 
    Ca = Ca_inf;
    E[1] = E[2] = 500.0*R_F*(T + 273.15)*log(3000.0/Ca);


  for (int ii = 0; ii < 7; ++ii)
  { 

    gate_inf[0][ii] = Xinf(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1]);
    if (ii == 4)
    {
      gate_inf[0][ii] = Xinf(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1],Ca/(Ca+3.0));
    }
    if (ii < 4)
    {
      gate_inf[1][ii] = Xinf(V,Xinf_par[ii][1][0],Xinf_par[ii][1][1]);
    }
    else
    {
      gate_inf[1][ii] = 1.0;
    }

    gate[0][ii] =  gate_inf[0][ii];
    gate[1][ii] =  gate_inf[1][ii];
    gnmh[ii] = 0.0;
  }




}

double neuron::Ca_dV( void )
{  
  /* This is the derivative of the Calcium Concentration with the following voltage dependency:
      [Ca] = - (0.94 * 10) * I_Ca(V)+0.05
      and: I_Ca = g_CaT * m_CaT^3 * h_CaT * (V-E_Ca) + g_CaS * m_CaS^3 * h_CaS * (V-E_Ca) 
      for each contribution the derivative can be written as:
      g_X (m_X^3 * h_X + [3 m_X^2 * h_X * dm_X/dV + m_X^3 * dh_X/dV] * (V-E_Ca))
      with dm_X/dV = Xinf_dV(V0,A_X,B_X) */

      // readout of the parameters and assigning of the dummy-variables

  // first for CaT
  float A_CaT_m = Xinf_par[1][0][0];
  float B_CaT_m = Xinf_par[1][0][1];
  float A_CaT_h = Xinf_par[1][1][0];
  float B_CaT_h = Xinf_par[1][1][1];
  double m_CaT = Xinf(V,A_CaT_m,B_CaT_m);
  double h_CaT = Xinf(V,A_CaT_h,B_CaT_h);
  double dm_CaT_dV = Xinf_dV(V,A_CaT_m,B_CaT_m);
  double dh_CaT_dV = Xinf_dV(V,A_CaT_h,B_CaT_h);

  // second for CaS
  float A_CaS_m = Xinf_par[2][0][0];
  float B_CaS_m = Xinf_par[2][0][1];
  float A_CaS_h = Xinf_par[2][1][0];
  float B_CaS_h = Xinf_par[2][1][1];
  double m_CaS = Xinf(V,A_CaS_m,B_CaS_m);
  double h_CaS = Xinf(V,A_CaS_h,B_CaS_h);
  double dm_CaS_dV = Xinf_dV(V,A_CaS_m,B_CaS_m);
  double dh_CaS_dV = Xinf_dV(V,A_CaS_h,B_CaS_h);

  // Calculation following the principle described above 
  double X = 0;
  X =  g[1] * ( (get_power(m_CaT,3))*h_CaT +( (3*get_power(m_CaT,2))*h_CaT*dm_CaT_dV + (get_power(m_CaT,3))*dh_CaT_dV ) * (V - E[1]) );
  X += g[2] * ( (get_power(m_CaS,3))*h_CaS +( (3*get_power(m_CaS,2))*h_CaS*dm_CaS_dV + (get_power(m_CaS,3))*dh_CaS_dV ) * (V - E[1]) );
  X *= - 0.94;

  return X;
}

double neuron::KCainf_dV( bool Ca_ = false )
{
    /* this function is the special version of function Xinf
        Here we look especially at the voltage dependency of the Constant C 
        So the style is like dC/dV * Xinf + dXinf/dV * C
        one step further we can split dC/dV * Xinf = dC/d[Ca] * d[Ca]/dV * Xinf 
        with d[Ca]/dV = Ca_dV */
    float A = Xinf_par[4][0][0];
    float B = Xinf_par[4][0][1];
    double X;

    if (Ca_){
      X =  (1./(Ca+3.)-Ca/(get_power((Ca+3.0),2))) *Xinf(V,A,B)*Ca_dV() ; 
    }
    else{
      X=  Xinf_dV(V,A,B)*Ca/(Ca+3.0);
    }
    return X;
}
 

double neuron::get_DIC(int ii_){
  return DIC[ii_];
}


void neuron::calculate_DIC( void )
{
  for (int jj = 0; jj < 3; ++jj)
  {
    DIC[jj] = 0;
  }
  double factor[3][2] = {{0,0},{0,0},{0,0}};
  double w_fs_m = 0;
  double w_fs_h = 0;
  double w_su_m = 0;
  double w_su_h = 0;

  for (int ii = 0; ii < 7; ++ii)
  {
    weight(ii);
    w_fs_m = weight_array[0];
    w_fs_h = weight_array[1];
    // the function returns [w_fs,w_su] with [m,h] for each
    w_su_m = weight_array[2];
    w_su_h = weight_array[3];

    // the first factor means w_fs
    factor[0][0] = w_fs_m;
    factor[0][1] = w_fs_h;
    // second factor means (w_su - w_fs) 
    factor[1][0] =  w_su_m-w_fs_m;
    factor[1][1] =  w_su_h - w_fs_h;
    // third factor means 1- w_su
    factor[2][0] = 1 - w_su_m;
    factor[2][1] = 1 - w_su_h;

    float Ca_vec[3];

    if (ii == 4)
    {
      weight(7); // calculate the weight for tau_Ca
      Ca_vec[0] =  0; // weight_array[0];
      Ca_vec[1] =  0; // weight_array[2]-weight_array[0];
      Ca_vec[2] =  1; // 1 - weight_array[2];

    }
    for (int jj = 0; jj < 3; ++jj)
    { 
      if (ii != 4){
        DIC[jj] -= factor[jj][0] * g[ii]*(p[ii]*get_power(gate[0][ii],(p[ii]-1)))*gate[1][ii] * (V - E[ii])* Xinf_dV(V,Xinf_par[ii][0][0],Xinf_par[ii][0][1])/Capacity;
      }
      else{
        DIC[jj] -= factor[jj][0] * g[ii]*(p[ii]*get_power(gate[0][ii],(p[ii]-1)))*gate[1][ii] * (V - E[ii])* KCainf_dV() /Capacity;
        DIC[jj] -= Ca_vec[jj] * g[ii]*(p[ii]*get_power(gate[0][ii],(p[ii]-1)))*gate[1][ii] * (V - E[ii])* KCainf_dV(true) /Capacity;
      }
      // skip the calculation for the derivative of h, if there is no contribution in h
      if (ii < 4){
        DIC[jj] -= factor[jj][1] * g[ii]*(get_power(gate[0][ii],p[ii])) * (V - E[ii])* Xinf_dV(V,Xinf_par[ii][1][0],Xinf_par[ii][1][1])/Capacity;
        
      }
    }

  }
}

void neuron::weight(int  ii_){

  // weight_array = [w_fs_m, w_fs_h,w_su_m, w_su_h] 
  // assigning the fast, slow, and ultraslow time constants
  double tau_f = tau[0][0]; // corresponds to NA - tau_m 
  double tau_s = tau[0][5]; // corresponds to Kd - tau_m
  double tau_u = tau[1][2]; // corresponds to CaS - tau_h
  // double tau_u = 110;//tau[0][6]; // corresponds to H - tau_m
  
  // built up the weights


  for (int jj = 0; jj < 2; ++jj){
    if (tau[jj][ii_] <= tau_f)
    {
          weight_array[jj+2*0] = 1;
          weight_array[jj+2*1] = 1;
    }
    else if (tau_f < tau[jj][ii_] && tau[jj][ii_] <= tau_s)
    {
          weight_array[jj+2*0] = ( log(tau_s)-log(tau[jj][ii_]) )/( log(tau_s)-log(tau_f) );
          weight_array[jj+2*1] = 1;
    }
    else if (tau_s < tau[jj][ii_] && tau[jj][ii_] <= tau_u)
    {
          weight_array[jj+2*0] = 0;
          weight_array[jj+2*1] = ( log(tau_u)-log(tau[jj][ii_]) )/( log(tau_u)-log(tau_s) );
    }
    else if  (tau[jj][ii_] > tau_u){
          weight_array[jj+2*0] = 0;
          weight_array[jj+2*1] = 0;
    }
  }
}



double neuron::get_V( void ){
	return V;
}

double neuron::get_Ca( void ){
	return Ca;
}
void neuron::set_V( double V_){
  V = V_;
}
double neuron::get_membrane_current( void ){
  I = g_L*(V - E_L);
  for (int ii = 0; ii < 7; ++ii){
  
    I += g[ii] * get_power(gate[0][ii], p[ii])*gate[1][ii]*(V-E[ii]);
  }

  return I;
}



int simulation(int first_cell_ , int second_cell_ , float V_start_1_, float V_start_2_,
          float g_from_cell_1_to_cell_2_, float g_from_cell_2_to_cell_1_ , 
          float V_syn_reverse_cell_1_, float V_syn_reverse_cell_2_,
          float synapse_V_half_from_0_to_1, float synapse_V_slope_from_0_to_1, float synapse_tau_from_0_to_1,  
          float synapse_V_half_from_1_to_0, float synapse_V_slope_from_1_to_0, float synapse_tau_from_1_to_0,  
          const char* Filename_V_, const char* Filename_sp_, bool supress_output_, bool supress_V_file_output, bool supress_Sp_file_output, bool check_the_properties,
          int max_con_to_change_cell_0[], float change_max_con_cell_0_to_[], int how_many_changes_in_cell_0_,
          int max_con_to_change_cell_1[], float change_max_con_cell_1_to_[], int how_many_changes_in_cell_1_,
          int simulation_time,
          float min_period_threshold, float max_period_threshold,
          float min_phase_difference_threshold, float max_phase_difference_threshold,
          float min_burst_exclusion_threshold
          ){

bool cell_one_spiked = false;
bool cell_two_spiked = false;



if (not supress_output_){
  cout << endl << "############################################################" << endl;
  cout << "Starting simulation with the following arguments:" << endl;
  cout << "first_cell: " << first_cell_ << "  second_cell: " << second_cell_ << endl;
  cout << "V_start = " << V_start_1_ << ", " << V_start_2_ << endl;
  cout << "g_synapse: From cell 1 to cell 2: " << g_from_cell_1_to_cell_2_ << "       And From cell 2 to cell 1: " << g_from_cell_2_to_cell_1_ << endl;
  cout << "Synapse: V_half = " << synapse_V_half_from_0_to_1 << "   V_slope = " << synapse_V_slope_from_0_to_1 << "    tau_syn = " << synapse_tau_from_0_to_1 << endl;
  cout << "Synapse: V_half = " << synapse_V_half_from_1_to_0 << "   V_slope = " << synapse_V_slope_from_1_to_0 << "    tau_syn = " << synapse_tau_from_1_to_0<< endl;
  cout << "saving to:" << endl;
  cout << Filename_V_ << endl;
  cout << Filename_sp_ << endl; 
}

neuron neuron_0( first_cell_, V_start_1_ , g_from_cell_2_to_cell_1_, V_syn_reverse_cell_1_, supress_output_, max_con_to_change_cell_0, change_max_con_cell_0_to_, how_many_changes_in_cell_0_);
cout << "to here";
neuron neuron_1( second_cell_, V_start_2_ , g_from_cell_1_to_cell_2_, V_syn_reverse_cell_2_, supress_output_, max_con_to_change_cell_1, change_max_con_cell_1_to_, how_many_changes_in_cell_1_);


float V_th = -30;
bool voltage_high[2] = {false,false};
float spike_time[2] = {0,0};



double current_time = 0;
double Svar[2];
double S_inf[2];


double V_half_from_0_to_1 = synapse_V_half_from_0_to_1;
double V_half_from_1_to_0 = synapse_V_half_from_1_to_0;
double V_slope_from_0_to_1 = synapse_V_slope_from_0_to_1;
double V_slope_from_1_to_0 = synapse_V_slope_from_1_to_0;
double tau_syn_from_0_to_1 = synapse_tau_from_0_to_1;
double tau_syn_from_1_to_0 = synapse_tau_from_1_to_0;




char V_file [999];
sprintf (V_file, "%s",Filename_V_);

FILE *Vfile;

if (!supress_V_file_output)
{
  Vfile = fopen(V_file,"w");
}


// Here we store the spike train in vectors for futher in-program evaluation
// But we can also generate Files for further steps


char Sp_file [999];
sprintf (Sp_file, "%s",Filename_sp_);
FILE *Spfile;
if (!supress_Sp_file_output)
{ 
  Spfile = fopen(Sp_file,"w");
  fclose(Spfile);
}

std::vector<float>Sp_vector_cell_0;
std::vector<float>Sp_vector_cell_1;
int number_of_spikes_cell_0 = 0;
int number_of_spikes_cell_1 = 0;

int total_time = simulation_time;
cout << "actually starting now";
while ( current_time < total_time){


    neuron_0.evolve_Ca(); neuron_1.evolve_Ca();
    neuron_0.calculate_time_constants(); neuron_1.calculate_time_constants();
    neuron_0.evolve_gates(); neuron_1.evolve_gates();
    // now the calculation for the synapses begins
    S_inf[1] = Xinf(neuron_0.get_V() , V_half_from_0_to_1, V_slope_from_0_to_1);
    S_inf[0] = Xinf(neuron_1.get_V() , V_half_from_1_to_0, V_slope_from_1_to_0);

    Svar[0] = Svar[0] + (1-exp(-dt/tau_syn_from_1_to_0))*(S_inf[0]-Svar[0]);
    Svar[1] = Svar[1] + (1-exp(-dt/tau_syn_from_0_to_1))*(S_inf[1]-Svar[1]);

    neuron_0.evolve_V(Svar[0]); 
    neuron_1.evolve_V(Svar[1]);

  if (!supress_V_file_output)
  {
  if ( 1000 > total_time-current_time && fmod(current_time,0.2)<dt)
  { 
    fprintf(Vfile,"%f %f %f\n", current_time,neuron_0.get_V(), neuron_1.get_V());
  }
  }
   if ( neuron_0.get_V() > V_th && current_time > 5000){
      if (!voltage_high[0] and (current_time - spike_time[0] > 5)){
          Sp_vector_cell_0.push_back(current_time);
          number_of_spikes_cell_0 ++;
          voltage_high[0] = true;
          spike_time[0] = current_time;
          cell_one_spiked = true;
      }
   }
   else{
      if (voltage_high[0]){
        voltage_high[0] = false;
      }
   }
  if ( neuron_1.get_V() > V_th && current_time > 5000){
      if (!voltage_high[1] and (current_time - spike_time[1] > 5)){
          Sp_vector_cell_1.push_back(current_time);
          number_of_spikes_cell_1 ++;
          voltage_high[1] = true;
          spike_time[1] = current_time;
          cell_two_spiked = true;
      }
   }
   else{
      if (voltage_high[1]){
        voltage_high[1] = false;
      }
   }
   

   current_time += dt;

}

if (!supress_V_file_output)
{
  fclose(Vfile);
}
if (!supress_Sp_file_output)
{ 
  Spfile = fopen(Sp_file,"a");
  for (int p = 0; p < number_of_spikes_cell_0; ++p)
  {
    fprintf(Spfile,"0 %f\n", Sp_vector_cell_0[p]);
  }
  for (int p = 0; p < number_of_spikes_cell_1; ++p)
  {
    fprintf(Spfile,"1 %f\n", Sp_vector_cell_1[p]);
  }
  fclose(Spfile);
}







if (cell_one_spiked == false || cell_two_spiked == false){
  if (!supress_Sp_file_output)
  {
    string str(Sp_file);
  remove(str.c_str());
  }
  
  if (!supress_V_file_output)
  {
    string str(V_file);
    remove(str.c_str());
  }

  Sp_vector_cell_0.clear();
  Sp_vector_cell_1.clear();

  if (cell_one_spiked == false)
  {
    return 1;
  }
  if (cell_two_spiked == false)
  {
    return 2;
  }
  


}
else{
  
bool fits_the_conditions;

if (check_the_properties){
fits_the_conditions = evaluate_spike_train(Sp_vector_cell_0, Sp_vector_cell_1, 
                                            number_of_spikes_cell_0 , number_of_spikes_cell_1 ,
                                            simulation_time, 
                                            min_period_threshold, max_period_threshold,
                                            min_phase_difference_threshold, max_phase_difference_threshold,
                                            min_burst_exclusion_threshold);

}
else{
  fits_the_conditions = true;
}
 Sp_vector_cell_0.clear();
 Sp_vector_cell_1.clear();



if (fits_the_conditions){

  return 0;

}
else{

  return 3;
}


}


}





int do_nothing(){
  int minus_one = -1;
  return minus_one;
}

float getThreshold(vector<float>& spike_times_, int length_of_spike_times)
{
  // Input is just a list of spike times 
  // output is the guessed threshold
  const int length_of_ISIs = length_of_spike_times - 1;
  float percentile_to_look_at = 0.9;
  const int percentile_entry = int(floor((1 - percentile_to_look_at) * length_of_ISIs));
  // This value is the place I need to look for the percentile value in a sorted array when 
  // looking from the end

  float ISIs[length_of_ISIs];
  float ISI_dummy;
  float minimum_ISI = spike_times_[1] - spike_times_[0];
  float percentile_ISI;
  float sp_thr;
  for (int ii = 0; ii < length_of_ISIs; ++ii)
  { 
    ISIs[ii] = 0;
  }
  for (int ii = 0; ii < length_of_ISIs; ++ii)
  { 
    ISI_dummy = spike_times_[ii+1] - spike_times_[ii];
    if (ISI_dummy < minimum_ISI){
      minimum_ISI = ISI_dummy;
    }
    
    for (int jj = 1; jj < percentile_entry + 2; ++jj){
      if (ISIs[length_of_ISIs - jj] < ISI_dummy){
        // If the new value is smaller than the actual value
        for (int kk = length_of_ISIs - (percentile_entry + 3); kk < (length_of_ISIs - jj); ++kk)
        {
          ISIs[kk] = ISIs[kk + 1];
        }
        ISIs[length_of_ISIs - jj] = ISI_dummy;
        
        break;
      }
      
    }
  }
  percentile_ISI = ISIs[length_of_ISIs - percentile_entry];
  
  sp_thr = (percentile_ISI + minimum_ISI)/2.;
  if (abs(sp_thr - ISIs[-1]) < 10 || abs(sp_thr- minimum_ISI) < 10)
  { 
    sp_thr = 0.99 * minimum_ISI;
  }


  return sp_thr;

}



bool evaluate_spike_train(vector<float>& spike_train_0, vector<float>& spike_train_1, 
              int length_of_spike_train_0 , int length_of_spike_train_1 ,
              int simulation_time_, 
              float min_period_threshold, float max_period_threshold,
              float min_phase_difference_threshold, float max_phase_difference_threshold,
              float min_burst_exclusion_threshold){
  // Here we ask for our conditions (phase differenc<e, burst exclusion and period)
  // True means that everything is ok
  

  /* READOUT THE SPIKE TRAINS */

      // get the Threshold for the burst detection
      float threshold_0 = getThreshold(spike_train_0,length_of_spike_train_0);


      float threshold_1 = getThreshold(spike_train_1,length_of_spike_train_1);

      // cout << endl << threshold_0 << " " << threshold_1 << endl;
      
  /* DETECT THE BURST FOR EACH INDIVIDUAL CELL */


        int nn = 0;


        // CELL 0
        /* leads to
        num_sp_0_cut
      first_spike_0_cut
      last_spike_0_cut 
      n_bursts_0 */
          float num_sp_0[length_of_spike_train_0];
          // This is the maximum amount possible for number of bursts
          float first_spike_0[length_of_spike_train_0];
          float last_spike_0[length_of_spike_train_0];
          int n_bursts_0 = 0;

          // initialize the first values
          num_sp_0[0] = 1;
          first_spike_0[0] = spike_train_0[0];

          while (length_of_spike_train_0 > nn + 1){

            if (spike_train_0[nn+1] - spike_train_0[nn] < threshold_0){
              // if spikes are close enough in time
              num_sp_0[n_bursts_0] ++;
            }
            else{
              last_spike_0[n_bursts_0] = spike_train_0[nn];
              n_bursts_0 ++;
              num_sp_0[n_bursts_0] = 1;
              first_spike_0[n_bursts_0] = spike_train_0[nn + 1];

            }
            nn ++;

          }
          last_spike_0[n_bursts_0] = spike_train_0[nn];

          float num_sp_0_cut[n_bursts_0];
          // This is the maximum amount possible for number of bursts
          float first_spike_0_cut[n_bursts_0];
          float last_spike_0_cut[n_bursts_0];

          for (int ii = 0; ii < n_bursts_0; ++ii)
          {
          num_sp_0_cut[ii] = num_sp_0[ii+1];
          first_spike_0_cut[ii] = first_spike_0[ii+1];
          last_spike_0_cut[ii] = last_spike_0[ii+1];
          }
        

        nn = 0;

        // CELL 1
        /* leads to
        num_sp_1_cut
      first_spike_1_cut
      last_spike_1_cut 
      n_bursts_1 */

          float num_sp_1[length_of_spike_train_1];
          // This is the maximum amount possible for number of bursts
          float first_spike_1[length_of_spike_train_1];
          float last_spike_1[length_of_spike_train_1];
          int n_bursts_1 = 0;

          // initialize the first values
          num_sp_1[0] = 1;
          first_spike_1[0] = spike_train_1[0];

          while (length_of_spike_train_1 > nn + 1){

            if (spike_train_1[nn+1] - spike_train_1[nn] < threshold_1){
              // if spikes are close enough in time
              num_sp_1[n_bursts_1] ++;
            }
            else{
              last_spike_1[n_bursts_1] = spike_train_1[nn];
              n_bursts_1 ++;
              num_sp_1[n_bursts_1] = 1;
              first_spike_1[n_bursts_1] = spike_train_1[nn + 1];

            }
            nn ++;

          }
          last_spike_1[n_bursts_1] = spike_train_1[nn];

          float num_sp_1_cut[n_bursts_1];
          // This is the maximum amount possible for number of bursts
          float first_spike_1_cut[n_bursts_1];
          float last_spike_1_cut[n_bursts_1];

          for (int ii = 0; ii < n_bursts_1; ++ii)
          {
          num_sp_1_cut[ii] = num_sp_1[ii+1];
          first_spike_1_cut[ii] = first_spike_1[ii+1];
          last_spike_1_cut[ii] = last_spike_1[ii+1];
          }

  /* DETECT THE BURSTS IN THE CIRCUIT */
  //  Now we got the arrays 
  //  first_spike_0, last_spike_0, num_sp_0
  //  with useful length n_ciruit_bursts_0
  //  AND
  //  first_spike_1, last_spike_1, num_sp_1
  //  with useful length n_ciruit_bursts_1
    /* 
       Here we want to capture the alternating effects of the circuit
       such that we want to call everything a burst, that happens between two bursts of the other cell
       In case there is an overlap of two bursts, this is counted as a new burst as well
    */
    float num_sp_A[max(n_bursts_1,n_bursts_0)];
    float num_sp_B[max(n_bursts_1,n_bursts_0)];
    float first_spike_A[max(n_bursts_1,n_bursts_0)];
    float first_spike_B[max(n_bursts_1,n_bursts_0)];
      float last_spike_A[max(n_bursts_1,n_bursts_0)];
    float last_spike_B[max(n_bursts_1,n_bursts_0)];
    float len_spike_train_A;
    float len_spike_train_B;

    /* PICTURE OF WHAT IS HAPPENING */

    /*
        This is considered to be            This is considered to be 
        the first burst of Trace A            the next burst of Trace A

            We call it                    We call it 
            "actual_A"                    "next_A"
        ##########################            ##########################
        

                        ##############################
                         This is considered to be the 
                             burst of trace B 
                            
                            We call this
                             "actual_B"

        Each Burst is characterized by a _beginning, and _ending, and a _num_sp
        except for the ones in the B-trace because it does not matter how many spikes we got
    */
    float actual_A_beginning;
    float actual_A_ending;
    float actual_A_num_sp;
    float next_A_beginning;
    float next_A_ending;
    float next_A_num_sp;
    float actual_B_beginning;
    float actual_B_ending;

    // Also we need an array to store the new burst characteristics in
    float first_spike[max(n_bursts_0,n_bursts_1)];
    float last_spike[max(n_bursts_0,n_bursts_1)];
    float num_sp[max(n_bursts_0,n_bursts_1)];
    int n_ciruit_bursts;

    int last_kk;


    /* important numbers after the assignment */
    int n_ciruit_bursts_0 = 0;
    int n_ciruit_bursts_1 = 0;
    bool found_next_burst = false;
    bool break_case;
    for (int omega = 0; omega < 2; ++omega)
    {
    // begin omega loop
      
      if (omega == 0)
      { 
        len_spike_train_A = n_bursts_0;
        len_spike_train_B = n_bursts_1;

        for (int ii = 0; ii < n_bursts_0; ++ii)
        {
          num_sp_A[ii] = num_sp_0_cut[ii];
          first_spike_A[ii] = first_spike_0_cut[ii];
          last_spike_A[ii] = last_spike_0_cut[ii];
        }
        for (int ii = 0; ii < n_bursts_1; ++ii)
        {
          num_sp_B[ii] = num_sp_1_cut[ii];
          first_spike_B[ii] = first_spike_1_cut[ii];
          last_spike_B[ii] = last_spike_1_cut[ii];
        }
      }
      if (omega == 1)
      { 
        len_spike_train_A = n_bursts_1;
        len_spike_train_B = n_bursts_0;
        for (int ii = 0; ii < n_bursts_0; ++ii)
        {
          num_sp_B[ii] = num_sp_0_cut[ii];
          first_spike_B[ii] = first_spike_0_cut[ii];
          last_spike_B[ii] = last_spike_0_cut[ii];
        }
        for (int ii = 0; ii < n_bursts_1; ++ii)
        {
          num_sp_A[ii] = num_sp_1_cut[ii];
          first_spike_A[ii] = first_spike_1_cut[ii];
          last_spike_A[ii] = last_spike_1_cut[ii];
        }
      }
      nn = 0;
      last_kk = 0;
      n_ciruit_bursts = 0;
      break_case = false;
      while (nn < len_spike_train_A)
      {
      // Start the loop through all the A-bursts


        actual_A_beginning = first_spike_A[nn];
        actual_A_ending = last_spike_A[nn];
        actual_A_num_sp = num_sp_A[nn];
        first_spike[n_ciruit_bursts] = first_spike_A[nn];
        last_spike[n_ciruit_bursts] = last_spike_A[nn];
        num_sp[n_ciruit_bursts] = num_sp_A[nn];

        for (int jj = nn + 1; jj < len_spike_train_A; ++jj)
        { 
        // begin loop NEXT
          /*
          In this loop we start from the next burst in trace A
          and will loop as long as we cross a beginning of a 
          burst in trace B 
          */

          
          next_A_beginning = first_spike_A[jj];
          next_A_ending = last_spike_A[jj];
          next_A_num_sp = num_sp_A[jj];

          for (int kk = last_kk; kk < len_spike_train_B; ++kk)
          {
          // begin loop B-trace
            actual_B_beginning = first_spike_B[kk];
            actual_B_ending = last_spike_B[kk];



            /* NOW TWO DIFFERENT CASES */


            /*  PICTURE

              1)


                ############              ##############
                    #############
                
                or 
              
                ############              ##############
                          #############
                
                or
                ############              ##############
                                #############

              1*)
                    ############
                  #############
                or

                          ############
                  #############

              2)

                ############              ##############
                                      #############
                
                or 
              
                ############              ##############
                                            #############

            
            */

            if (actual_B_beginning < next_A_beginning )
            {
              if (actual_A_beginning < actual_B_beginning)
              { 
                // This is case 1
                last_kk = kk;
                found_next_burst = true;
                break;

              }
              else{
                // This is case 1*
                /* 
                If from now on we only see A - bursts 
                this basically means that there are no more B bursts coming after this A
                burst and thus this A burst should b neglected */
                if (kk == len_spike_train_B-1){
                  found_next_burst = false;
                  break;
                }

              }
            }
            else{ 
                
                // This is case 2
                last_kk = kk;
                found_next_burst = false;

                if (jj == len_spike_train_A - 1){
                  if (last_kk < len_spike_train_B - 1){
                    last_spike[n_ciruit_bursts] = next_A_ending;
                    num_sp[n_ciruit_bursts] += next_A_num_sp;
                    n_ciruit_bursts ++;
                    break_case = true;
                  }
                }



                break;
              
              
            }



          // end loop B-trace
          }


          if (found_next_burst)
          {
            nn = jj - 1;
            found_next_burst = false;
            n_ciruit_bursts += 1;
            break;
          }
          else{
            last_spike[n_ciruit_bursts] = next_A_ending;
            num_sp[n_ciruit_bursts] += next_A_num_sp;
          }

          
      
        // end loop NEXT
        }

       if (break_case){
        break;
       }

      nn += 1;


      // End the loop through all the A-bursts
      }

    // Finally we rewrite the beginning and end of the burst
    /* 

    !!!

      Attention here
      The arrays can be longer now than they are filled 
      with useful stuff 
      The important numbers will be stored in

      n_bursts_0 and n_bursts_1 

    !!!

    */

    if (omega == 0)
    {
      n_ciruit_bursts_0 = n_ciruit_bursts - 1;
      for (int ii = 0; ii < n_ciruit_bursts_0; ++ii)
      {
        first_spike_0[ii] = first_spike[ii+1];
        last_spike_0[ii] = last_spike[ii+1];
        num_sp_0[ii] = num_sp[ii+1];
      }
    }
    if (omega == 1)
    {
      n_ciruit_bursts_1 = n_ciruit_bursts - 1;
      for (int ii = 0; ii < n_ciruit_bursts_1; ++ii)
      {
        first_spike_1[ii] = first_spike[ii+1];
        last_spike_1[ii] = last_spike[ii+1];
        num_sp_1[ii] = num_sp[ii+1];
      }
    }

    // end omega loop
    }


    /* 

    Now we got the arrays 

    first_spike_0, last_spike_0, num_sp_0
    with useful length n_ciruit_bursts_0
    AND
    first_spike_1, last_spike_1, num_sp_1
    with useful length n_ciruit_bursts_1
    */
  

  /* controling output */
    /*
    cout << " 0 " << endl;
    for (int i = 0; i < n_ciruit_bursts_0; ++i)
    {
      cout << first_spike_0[i] << " " << last_spike_0[i] << " " << num_sp_0[i] << endl;
    }
    cout << n_ciruit_bursts_0;
    cout << " 0 " << endl;
    for (int i = 0; i < n_ciruit_bursts_1; ++i)
    {
      cout << first_spike_1[i] << " " << last_spike_1[i] << " " << num_sp_1[i] << endl;
    }
    */

  /* START WITH THE ANALYSIS NOW */

    /* BURST EXCLUSION METRIC */
    /*  this leads to chi
      the burst exclusion measure
    */


      // If one cell is only spiking, there is no 
      // further interest in this cell-pair
      int cell_0_is_only_spiking = 0;
      int cell_1_is_only_spiking = 0;
      for (int ii = 0; ii < n_ciruit_bursts_0; ++ii)
      {
        if (num_sp_0[ii] > 1){
          cell_0_is_only_spiking += 1;
          if (cell_0_is_only_spiking > 5){
            break;
          }
        }

      }
      for (int ii = 0; ii < n_ciruit_bursts_1; ++ii)
      {
        if (num_sp_1[ii] > 1){
          cell_1_is_only_spiking += 1;
          if (cell_1_is_only_spiking > 5){
            break;
          }
        }

      }

      if ( !(cell_0_is_only_spiking > 5 and cell_1_is_only_spiking > 5 )){
        return false;
      }

      /*
      # first case:
      #   s0#####e0 
      #            s1######e1

      #second case:
      #       s0####e0
      #   s1######e1
      */

      int last_jj = 0;
      float Onet = 0;


      float s0, e0;
      float s1, e1;


      for (int ii = 0; ii < n_ciruit_bursts_0; ++ii)
      {
      // Start the loop through trace 0
        s0 = first_spike_0[ii];
        e0 = last_spike_0[ii];

        if(s0 == e1){
          continue;
        }
        for (int jj = last_jj; jj < n_ciruit_bursts; ++jj)
        { 
        // Start the loop through trace 1
          s1 = first_spike_1[jj];
          e1 = last_spike_1[jj];

          
          if (s1 > e0){
            /*
              # here the first loop needs to catch the second one
              # so no 'last_jj'-update because we still want to use
              # the actual burst 1 for the next comparison 
            */
            break;    
          }
          else{
            last_jj = jj;
            if (e1 < s0)  
            {
              /* 
              # now the second loop needs to catch the first, so we just 
              # keep it running 
              */
              continue;
            }
            else{
              /*
              # now e1 >= s0 and e0 >= s1
              # possible cases:
              
              1: 
                  s0########e0
                s1########e1
              2:
                  s0#######e0
                    s1########e1
              3:    
                  s0########e0
                     s1###e1
              4:
                  s0########e0
                s1################e1
              */

              float burst_0 = abs(e0-s0);
              float burst_1 = abs(e1-s1);
              if ( (s0 > s1 && e0 < e1) || (s0 <= s1 && e0 > e1) )
              // cases 3 and 4
              { 
                Onet += min(burst_0,burst_1);
              }
              else{
              // cases 1 and 2
                if (s0 > s1){
                  Onet += burst_0 - (e0-e1);
                }
                else{
                  Onet += burst_1 - (e1-e0);
                }
              }
            }
          }
          if (e1 < e0)
          {
            /*
            # here it could be, that there is a second burst of line 1, 
            # which also intersects with burst 0 
            # so we just keep circulating in loop 1
            */
            continue;
          }
          if (e1 > e0)
          {
            /*
            # now here burst 1 ends after burst 0, so 
            # maybe there is another
            */
            break;
          }
        // Stop the loop through trace 1
        }
      // Stop the loop through trace 0
      }


      float sum_time_0 = 0;
      float sum_time_1 = 0;

      for (int ii = 0; ii < n_ciruit_bursts_0; ++ii)
      {
        sum_time_0 += (last_spike_0[ii] - first_spike_0[ii]);
      }

      for (int ii = 0; ii < n_ciruit_bursts_1; ++ii)
      {
        sum_time_1 += (last_spike_1[ii] - first_spike_1[ii]);
      }

      float T_trial = simulation_time_;


      float chi, Omin, Omax, Orand;

      if (sum_time_0 + sum_time_1 > T_trial)
      { 
        // overlap time if random
        Orand = min(sum_time_0,sum_time_1) - 0.5*(T_trial-max(sum_time_0,sum_time_1));
        // minimum possible overlap time
          Omin = sum_time_0 + sum_time_1 - T_trial;
      }
      else{
        // overlap time if random
        Orand = pow (min(sum_time_0,sum_time_1),2)/(2*(T_trial-max(sum_time_0,sum_time_1)));
        // minimum possible overlap time
          Omin = 0;
      }
      // maximum possible overlap
      Omax = min(sum_time_0,sum_time_1);
      if (Onet > Orand)
      {
        chi = (Orand-Onet)/(Omax-Orand);
      }
      else{
        chi = (Orand-Onet)/(Orand-Omin);
      }

    /* BURST EXCLUSION METRIC */

    /* MEAN PHASE DIFFERENCE AND MEAN PERIOD */
    /* This leads to 
      phase_difference and
      period
    */

      last_jj = 0;
      float phase_difference = 0;
      float period;

      float phase_difference_sum = 0;
      float period_sum = 0;

      float num_of_phase_differences = 0;
      float num_of_periods = 0;

      for (int ii = 0; ii < n_ciruit_bursts_0-1; ++ii)
      { 

        s0 = first_spike_0[ii];
        period_sum += first_spike_0[ii + 1] - s0;
        num_of_periods ++;
        for (int jj = last_jj; jj < n_ciruit_bursts_1; ++jj)
        {
          s1 = first_spike_1[jj];
          if (s1 < s0){
            continue;
          }
          else{
            last_jj = jj;
            phase_difference_sum += ((s1 - s0)/(first_spike_0[ii+1] - s0));
            num_of_phase_differences ++;
            break;

          }
        }
      }
      phase_difference = phase_difference_sum/num_of_phase_differences;
      period = period_sum/num_of_periods;

    /* MEAN PHASE DIFFERENCE AND MEAN PERIOD */


  /* EVALUATE WHETHER THIS FITS THE EXPECTATIONS */

    /* Now we need to qualify whether we want to keep this circuit 
    depending on the measure we just defined */
    if ( (period > min_period_threshold) && (period < max_period_threshold)  )
    { 
        if ((phase_difference > min_phase_difference_threshold) && (phase_difference < max_phase_difference_threshold))
        {
        if (chi > min_burst_exclusion_threshold)
        { 
          return true;
        }
      }
    }
    return false;



}







#endif

 


