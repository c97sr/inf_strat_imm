Last updated 07 Feb 2015

Task1: A clean version of serological model using jave component

Hong kong demographic: http://michellebian3.files.wordpress.com/2011/10/by-age.png

Files structure
- main
 - calculate_R0
 - extract_antibody_titres  (pre-processing)
 - immune_boosting          (viz)
 - likelihood_estimation    (likelihood_cal)
 - mcmc                      
 - mplot
 - peaktime
 - produce_samples
 



	- calculate_R0
		calR0_byT  	  The function calls cal_NGM_bybeta2
		calR0_bypar 	  The function calls cal_NGM_bybetaByT

	- extract_antibody_titres
		main_producetitres_hk

	- immune_boosting
		plot_boosting
		plot_susc_2pars
        	plot_susc_2pars_S2 (paper)
        	plot_susc_age (among different age grous)       

	- likelihood_estimation
		main_pointestimatelikelihood_hk_java
		main_findmin_sera_age_java


	- mcmc
		calDIC
	        main_MCMC.m
        

- lib
	- chart
	- model
		- llh
			calculateLogLikelihood
			estimatelikelihood
			get_titres_prob (<-estimatelikelihood.m)
		- rt
			cal_NGM_bybetaByT(beta,Ab,pars)  Calclate NGM by whole sampling periods
			cal_NGM_bybeta2(beta,Ab,pars)	 Calclate NGM by sampling time T1 and T2
			cal_NGM
			getMOF_byS (called by cal_NGM.m)
		- susc
			getSusc
			make_h
			make_g

		InitParameters
		getImmBoost
		getMOF
		getTotalInfect.m (<-simulate)
		make_ics
		make_ics_fromtitres
		make_ics_fromtitres_byage
		make_ics_naive
		odef_islmod
		odef_islmodjava
		setParameters
		simulate
		_genBoostingRate.m (unclear)
		_getSerology (unclear)
		_my_estimatelikelihood (unclear)
  
	- optim
		getLLH
		getNegLLH
		getNegLLHAge
		getNegLLHAgejava*
		getNegLLHjava


	- sys
		runMCMC

Flow:
* main scripts:
1. extract_antibody_titres/main_producetitres_hk:		        extract and save antibody titres
2. likelihood_estimation/main_pointestimatelikelihood_hk_java: 	read the samples and calculate the likelihood surface
3. likelihood_estimation/main_findmin_sera_age_java:            estimate ML using java component
4. mcmc/main_MCMC:                                              estimate posterior using MCMC


main_producetitres_hk
-> extract_titres (extract_titres_table_paired.m)
-> getColumnID
-> save 'h1n1_titres.mat'

main_pointestimatelikelihood_hk_java
-> simulate (lib/model)	
-> odef_islmodjava (lib/model)
-> estimatelikelihood (lib/model)

main_findmin_sera_age_java
-> fmincon 
-> getNegLLHAgejava
-> matlabjava.jar
-> save('maxllh.mat')

main_MCMC
-> getNegLLHAgejava
-> matlabhavjava.jar
-> save('mcmc_6pars.mat')

% Figures
1. Comparison of titre model fit
2. Disease and serological dynamics of the titre model simulation

Figure name/ file name


* plot scripts:
main_plot_seroloy
plot_dynamics
plot_dynamics_byage: an alternative function which plot same activities as plot_dynamics but with different age groups.
plot_dynamics_frompost: plot dynamics using posterior distribution
plot_dynamics_immesbyage: an alternative immune boosting mechanism
plot_lhsurf: plot likelihood surface


* plot immune boosting:
plot_boosting
plot_boosting_immes: boosting by immune escpape mechanism.
plot_susc_2pars: plot immune protection


* output
project names
1. herdimmunity
   output: titres.mat, likelihood_output.mat, figures
2. herdimmunity_hk
   output: titres.mat, likelihood_output.mat, figures


* Others:
1. how to do model selection
check http://en.wikipedia.org/wiki/Deviance_information_criterion
store log likelihood for MCMC

2. build a java component to run disease model
metropolis_multipars_main_java
build class and export jre.
add java path
http://www.mathworks.co.uk/help/matlab/matlab_external/bringing-java-classes-and-methods-into-matlab-workspace.html#f111131
http://stackoverflow.com/questions/9520503/calling-java-from-matlab
PosteriorSamples stores all the samples from posterior distribution



