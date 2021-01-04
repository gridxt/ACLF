#pragma once
struct Ybus_Struct
{
	
	/* Configuration Parameters */
	int CoreEngine; // To select computation Engine: 0: MKL-DENSE, 1: MKL-SP-CSR	
	int Num_Threads; // Number of threads to be used 
	int OMP_Enable; // Enable or Disable OMP: 0: Disable 1: Enable
	int Enable_Memory_Stat; // 0: Do not gather memory usage statistics by MKL, 1: Gather memory usage statisyics by MKL
		
	// Store memory usage statistics in the following variables
	double Peak_Mem_Used; // peak memory used by MKL in KB	

	/* These variables hold Yf and Yt values in dense format which are used in PF_Solution as well */
	MKL_Complex16 *Yt_Dense;
	MKL_Complex16 *Yf_Dense;

	/* These variables hold Yf and Yt values in sparse format which are used in PF_Solution as well */
	MKL_Complex16 *Yf_Val;
	MKL_INT *Yf_IA;
	MKL_INT *Yf_JA;
	MKL_Complex16 *Yt_Val;
	MKL_INT *Yt_IA;
	MKL_INT *Yt_JA;

	/* These variables hold Cf and Ct values in dense format which are used in PF_Solution as well */
	MKL_Complex16 *CtPrime_Dense;
	MKL_Complex16 *CfPrime_Dense;

	/* These variables hold Cf and Ct values in sparse format which are used in PF_Solution as well */
	sparse_matrix_t Cf_CSR_Handle;
	MKL_Complex16 *Cf_Val;
	MKL_INT *Cf_IA;
	MKL_INT *Cf_JA;
	sparse_matrix_t Ct_CSR_Handle;
	MKL_Complex16 *Ct_Val;
	MKL_INT *Ct_IA;
	MKL_INT *Ct_JA;
	/* These variables hold Cf and Ct values in sparse format which are used in PF_Solution as well */
	MKL_Complex16 *Tap_Dense;
	sparse_matrix_t Tap_CSR_Handle;
	MKL_Complex16 *Tap_Val;
	MKL_INT *Tap_IA;
	MKL_INT *Tap_JA;
		
	/* Series admittance */
	MKL_Complex16 *Ysc;

	/* These variables are used in PFsolution only */
	MKL_Complex16 *A_Dense;	
	sparse_matrix_t A_CSR_Handle;	
	MKL_INT *A_IA;
	MKL_INT *A_JA;
	MKL_INT	*pointerE_A; //useless but needed!
	MKL_Complex16 *A_Val;

	/* This variable holds the final Ybus in dense format */
	MKL_Complex16 *Ybus_Dense;

	/* These variable hold the final Ybus in sparse format  */
	sparse_matrix_t Ybus_CSR_Handle; // Ybus handle in CSR format	
	MKL_INT *Ybus_IA;
	MKL_INT	*Ybus_JA;
	MKL_INT	*pointerE_Ybus; //useless but needed!
	MKL_Complex16  *Ybus_Values;
};







struct PFSol_Struct
{	
	//Timer to keep the execution time
	double *Timer;

	// Convert the volatge angles to degree as well
	double *Vm;
	double *Vang_deg;

	// Generation at each bus
	double *Pg;
	double *Qg;
	
	// Branch power flows
	double *PF;
	double *QF;
	double *PT;
	double *QT;

	// Branch losses
	MKL_Complex16 *Loss;
	double *FCHG;
	double *TCHG;

	//Partial derivatives
	MKL_Complex16 *dloss_dV_a_Dense;	
	sparse_matrix_t dloss_dV_a_CSR_Handle;
	MKL_INT *dloss_dV_a_IA;
	MKL_INT *dloss_dV_a_JA;
	MKL_INT *dloss_dV_a_pointerE;
	MKL_Complex16 *dloss_dV_a_val;

	MKL_Complex16 *dloss_dV_m_Dense;
	sparse_matrix_t dloss_dV_m_CSR_Handle;
	MKL_INT *dloss_dV_m_IA;
	MKL_INT *dloss_dV_m_JA;
	MKL_INT *dloss_dV_m_pointerE;
	MKL_Complex16 *dloss_dV_m_val;

	MKL_Complex16 *dchg_dVm_f_Dense;
	sparse_matrix_t dchg_dVm_f_CSR_Handle;
	MKL_INT *dchg_dVm_f_IA;
	MKL_INT *dchg_dVm_f_JA;
	MKL_INT *dchg_dVm_f_pointerE;
	MKL_Complex16 *dchg_dVm_f_val;
	
	MKL_Complex16 *dchg_dVm_t_Dense;
	sparse_matrix_t dchg_dVm_t_CSR_Handle;
	MKL_INT *dchg_dVm_t_IA;
	MKL_INT *dchg_dVm_t_JA;
	MKL_INT *dchg_dVm_t_pointerE;
	MKL_Complex16 *dchg_dVm_t_val;

	/* *** Mis. information *** */
	// Bus number of maximum and minimum voltage magn. and angle
	MKL_INT Vmin_mag_bus;
	MKL_INT Vmax_mag_bus;
	MKL_INT Vmin_ang_bus;
	MKL_INT Vmax_ang_bus;
	
	// Branch loss information	
	double Q_loss_total;
	double P_loss_total;
	MKL_INT P_loss_max_line_ind;
	MKL_INT Q_loss_max_line_ind;

	// Generation capacity/actual
	double Pg_total;
	double Qg_total_min;
	double Qg_total_max;
	double Pg_total_online;
	double Qg_total_min_online;
	double Qg_total_max_online;	
	double Pg_actual;
	double Qg_actual;

	// Load
	double Pd_total;
	double Qd_total;

	// Shunt injection
	double Shunt_P_inj;
	double Shunt_Q_inj;
	MKL_INT Shunt_number;

	// Branch charging
	double Branch_Q_inj;

	// Number of buses
	MKL_INT Bus_number;

	// Number of In- and out-of-service branches
	MKL_INT Branch_in_number;
	MKL_INT Branch_out_number;

	// Number of online and offline generators
	MKL_INT Gen_total_number;
	MKL_INT Gen_online_number;
	MKL_INT Gen_offline_number;

	// Number of transformers
	MKL_INT Transformers_number;

	// Number of loads
	MKL_INT Loads_number;

	// Status
	int status;

	// Memory usage
	double Peak_Mem_Usage;

};

