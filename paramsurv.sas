/*====================================================================================================================================================================
Parmsurv: a SAS Macro for Flexible Parametric Survival Analysis with Long-Term Predictions

Copyright © 2022 Merck Sharp & Dohme LLC., a subsidiary of Merck & Co., Inc., 126 East Lincoln Ave., P.O. Box 2000 Rahway, New Jersey 07065, USA. All rights reserved. 
Authors: 
	Han Fu, The Ohio State University, Columbus, OH, USA
	Shahrul Mt-Isa, MSD, Zurich, Switzerland
	Richard Baumgartner, Merck & Co., Inc., Rahway, NJ, USA
	Harris Kampouris, MSD, London, United Kingdom
	William Malbecq, MSD Europe INC., Brussels, Belgium

Parameters: 
	Dataset-related: data
	Response-related: t1, t2, censor, censval
	Covariate-related: covars, anc, class_cov, refgrp 
	Distribution: dist
	Optimization-related: optim_method, init, lower, upper
	Inference-related: alpha, robust
	Custom distribution: density, survival, hazard, custom_prep, nana, location, param_anc, param_anc_transf, log_transf_index, log_density, log_survival
	Prediction-related: pred, pred_max_time, pred_plot_cl
	Others: weight, strata, noprint
	(Please see documentation for detailed description of the parameters)

Statement: You may contact Han Fu at fu.607@osu.edu if you have questions or feedbacks/suggestions, but technical support or maintenance is not guaranteed.

Disclaimers: 
Han Fu was employed as an intern at Merck Sharp & Dohme LLC., a subsidiary of Merck & Co., Inc. during the study conduct. 
Shahrul Mt-Isa, Richard Baumgartner and William Malbecq are employees of Merck Sharp & Dohme LLC., a subsidiary of Merck & Co., Inc.

Users may use the codes provided freely at their own risks. Technical and maintenance support are not available and cannot be guaranteed. 
Although the codes have been validated, Authors and their affiliated organisations do not take responsibility of the results generated from these codes.

Additional details are available at https://arxiv.org/abs/2106.14109.
====================================================================================================================================================================*/


options mprint symbolgen mlogic;


%macro param_prep()/minoperator;
/*Function: connect optimization parameter theta to parameters needed in likelihood, used in each step of 
NLP optimization*/	
	theta_index = 1;
	param_value = j(ncol(param_vec),1); *param_vec = {"beta" "sigma" "lambda"} in GG;
	do param_id=1 to ncol(param_vec);
		index = ListGetItem(cov_index, param_vec[param_id]);
		if ncol(index)=1 then do; *no covariates associated with this parameter;
			param_value[param_id]=theta[theta_index];
			theta_index = theta_index + 1;
		end;
		else do; *there are covariates associated with this parameter;
	    	coef = theta[theta_index:(theta_index+ncol(index)-1)];
	    	param_value[param_id] = x_i[1,index]*coef;
			theta_index = theta_index + ncol(index);
		end;
	end;
	%if &dist ne %then %do; *built-in distributions;
		beta = param_value[1];
		%if %upcase(&dist.)=GENGAMMA %then %do;
			log_sigma = param_value[2];
			sigma = exp(log_sigma);
			lambda = param_value[3];
	        %if lambda>0 %then
		        %let survival=&survival1;
	        %else %let survival=&survival2;
	    %end; 
		%if %upcase(&dist.)=GENF %then %do;
			log_sigma = param_value[2];
			sigma = exp(log_sigma);
			q = param_value[3];
			log_p = param_value[4];
			p = exp(log_p);
		    d = sqrt(q*q+2*p);
		    m1 = 2/(q*q + 2*p + q*d);
		    m2 = 2/(q*q + 2*p - q*d);
		%end;
		%if %upcase(&dist.)=GENF_ORIG %then %do;
			log_sigma = param_value[2];
			sigma = exp(log_sigma);
		    log_m1 = param_value[3];
			m1 = exp(log_m1);
		    log_m2 = param_value[4];
			m2 = exp(log_m2);
		    d = sqrt(1/m1+1/m2);
		%end;
		%if %substr(%upcase(&dist.),1,3) = EXP %then %do;
			gamma = exp(-beta);
		%end;
		%if %upcase(&dist.) = GOMPERTZ %then %do;
			gamma = exp(-beta);
			shape = param_value[2];
		%end;
		%if %upcase(&dist.) in WEIBULL GAMMA LNORM LLOGIS %then %do;
			log_sigma = param_value[2];
			sigma = exp(log_sigma);
		%end;
	%end;
	%else %do; *custom distributions;
		%str(&exec);
		%str(&custom_prep);
	%end;
%mend param_prep;


%macro class_cov_prep (dataset=,where_st=)/minoperator;
/*Function: Prepare covariate matrix for calculation, including create dummy variables for classification covariates,
create a list that contains covariate labels and indices*/
	x_table = TableCreateFromDataSet("Work", "&dataset", "&where_st");
	nobs = nrow(x_table);
	x_num = {};
	refgrp_str = "&refgrp"; *reference group for classificaiton covariates;
	char_id = 0;
	ncov_com = ncol(covars_combined); *covars_combined is a vector of unique combined covariates from &covars and &anc;
	%if &class_cov ne %then %do;
		%let id = 1;
		class_cov_vec = {};
		%do %while (%scan(&class_cov, &id) ne );
			%let class_cov_s = %scan(&class_cov, &id);
			class_cov_vec = class_cov_vec || "&class_cov_s";
			%let id = %eval(&id+1);
		%end;
	%end;
	cov_name_dic = ListCreate(ncov_com); * list that stores covariate lables and indices for each covariate;
	call ListSetName(cov_name_dic,1:ncov_com,covars_combined);
	index = 1;
	do cov_id=1 to ncov_com; *for each covaraite;
		covar = covars_combined[cov_id]; *covariate name ;
		xc = TableGetVarData(x_table, covar); *covariate column;
		if %eval(&class_cov ne) then fl = any(class_cov_vec=covar); *flag: whether covariate specified in &class_cov;
		else fl=0;
		if type(xc)='N' & fl=0 then do; *continuous covariate;
			x_num = x_num || xc; 
			cov_dic_sub = ListCreate({"covarLabel" "covarIndex"});
			call ListSetItem(cov_dic_sub,"covarLabel",covar);
			call ListSetItem(cov_dic_sub,"covarIndex",index);
			call ListSetItem(cov_name_dic,covar,cov_dic_sub);
			index = index + 1;
		end;
		else do; *classification covariate;
			if type(xc)='N' then xc = putn(xc,"BEST6.");
			char_id = char_id + 1;
			category = unique(xc);
			call sort(category);
			ncat = ncol(category);
			if %eval(%nrstr(&refgrp)=) then refgrp = category[1];
			else refgrp = scan(refgrp_str,char_id,',');
			catded = setdif(upcase(strip(category)), upcase(strip(refgrp))); *categories without reference group;
			dummy = j(nobs, ncat-1, 0); *dummy variables;
			cov_name = {};
			do duco = 1 to (ncat-1);
				dummy[,duco] = (upcase(strip(xc))=strip(catded[duco]));
				cov_name = cov_name // catx("_",covar,lowcase(catded[duco]));
			end;
			cov_dic_sub = ListCreate({"covarLabel" "covarIndex"});
			call ListSetItem(cov_dic_sub,"covarLabel",cov_name);
			call ListSetItem(cov_dic_sub,"covarIndex",do(index,index+ncat-2,1));
			call ListSetItem(cov_name_dic,covar,cov_dic_sub);
			index = index + ncat - 1;
			x_num = x_num || dummy;
		end;
	end;
%mend class_cov_prep;

%macro sublist(param=,covars=);
/*Function: create sublist of param_list for each parameter and add the sublist to param_list*/
		covars_vec = {};
		%let cov_id = 1;
		%do %while (%scan(&covars, &cov_id) ne );
			%let covar = %scan(&covars, &cov_id);
			covars_vec = covars_vec || "&covar";
			%let cov_id = %eval(&cov_id+1);
		%end;
		if ncol(covars_vec)>0 then covars_com_vec = covars_com_vec || covars_vec;
		call ListSetItem(param_list, &param,covars_vec);
%mend sublist;



%macro paramsurv (data=,t1=,t2=,censor=,censval=0,covars=,anc=,class_cov=,refgrp=,weight=,strata=,
                    optim_method=nlpnra,dist=,init=,alpha=0.05,robust=no,nlp_print=0,log_result=no,
                    lower = {. . . . . . . . . .},upper = {. . . . . . . . . .},
                    density=,survival=,hazard=,custom_prep=,location=beta,param_anc=,log_transf_param=,log_density=,log_survival=,
                    pred=,pred_max_time=,pred_plot_cl=yes,res_print=yes)/minoperator;

	%let log0 = %str(10**(-15));


    %****** Define survival and density functions ******;

	%if &dist ne %then %do;
	    %if %upcase(&dist.)=GENGAMMA %then %do;
	        %let survival1=
	            %str(1-max(min(cdf('gamma',time**(lambda/sigma),1/(lambda*lambda),lambda*lambda*exp(beta*lambda/sigma)),1-&log0),&log0));
	        %let survival2=
	            %str(max(min(cdf('gamma',time**(lambda/sigma),1/(lambda*lambda),lambda*lambda*exp(beta*lambda/sigma)),1-&log0),&log0)));
	        %let density=%str(max(pdf('gamma',time**(lambda/sigma),1/(lambda*lambda),lambda*lambda*exp(beta*lambda/sigma)) * time**(lambda/sigma) * abs(lambda) / sigma,&log0));
			%let log_density = %str(log(max(pdf('gamma',time**(lambda/sigma),1/(lambda*lambda),lambda*lambda*exp(beta*lambda/sigma)),&log0)) + lambda/sigma*log(time) + log(max(abs(lambda),&log0)) - log_sigma);
	        %let param_anc = sigma lambda;
			%let log_transf_param = sigma;
	    %end;
	    %else %if %substr(%upcase(&dist.),1,4)=GENF %then %do;
	        %let survival=%str(cdf('beta',m2/(m2+m1*exp(d*(log(time)-beta)/sigma)),m2,m1));
	        %let density=%str(d*exp(m1*d*(log(time)-beta)/sigma)*(m1/m2)**m1/sigma/beta(m1,m2)/(1+m1*exp(d*(log(time)-beta)/sigma)/m2)**(m1+m2));
			%let log_density = %str(log(max(d,&log0)) + m1*d*(log(time)-beta)/sigma + m1*log(max(m1/m2,&log0)) - log_sigma - log(max(beta(m1,m2),&log0)) - (m1+m2)*log(max(1+m1*exp(d*(log(time)-beta)/sigma)/m2,&log0)));
			%if %upcase(&dist.)=GENF %then %do;
	        	%let param_anc = sigma q p;
				%let log_transf_param = sigma p;
			%end;
	        %else %do;
	        	%let param_anc = sigma m1 m2;
				%let log_transf_param = sigma m1 m2;
			%end;
	    %end; 
		%else %if %substr(%upcase(&dist.),1,3)=EXP %then %do;
			%let survival=%str(exp(-gamma*time));
			%let log_survival = %str(-gamma*time);
			%let density=%str(gamma*exp(-gamma*time));
			%let log_density = %str(-beta-gamma*time);
		%end;
		%else %if %upcase(&dist.)=GOMPERTZ %then %do;
			%let survival=%str(exp(-gamma/shape*(exp(shape*time)-1)));
			%let log_survival = %str(-gamma/shape*(exp(shape*time)-1));
			%let density=%str(gamma*exp(shape*time)*exp(-gamma/shape*(exp(shape*time)-1)));
			%let log_density = %str(-beta+shape*time-gamma/shape*(exp(shape*time)-1));
			%let param_anc=shape;
		%end;
		%if %upcase(&dist.)=WEIBULL %then %do;
	        %let survival=%str(exp(-exp((log(time)-beta)/sigma)));
	        %let log_survival=%str(-exp((log(time)-beta)/sigma));
	        %let density=%str(max(pdf('expo',time**(1/sigma),exp(beta/sigma)) * time**(1/sigma) / sigma,&log0));
			%let log_density = %str(log(max(pdf('expo',time**(1/sigma),exp(beta/sigma)),&log0)) + log(time)/sigma - log_sigma);
		%end;
		%if %upcase(&dist.)=GAMMA %then %do;
	        %let survival=%str(1-max(min(cdf('gamma',time,1/(sigma*sigma),sigma*sigma*exp(beta)),1-&log0),&log0));
	        %let density=%str(max(pdf('gamma',time,1/(sigma*sigma),sigma*sigma*exp(beta)) * time,&log0));
			%let log_density = %str(log(max(pdf('gamma',time,1/(sigma*sigma),sigma*sigma*exp(beta)),&log0)) + log(time));
		%end;
		%if %upcase(&dist.)=LNORM %then %do;
	        %let survival=%str(1-max(min(cdf('logn',time,beta,sigma),1-&log0),&log0));
	        %let density=%str(max(pdf('logn',time,beta,sigma),&log0));
		%end;
		%if %upcase(&dist.)=LLOGIS %then %do;
	        %let survival=%str(1/(1+exp(sqrt(2)*(log(time)-beta)/sigma)));
	        %let density=%str(sqrt(2)*exp(sqrt(2)*(log(time)-beta)/sigma)/sigma/(1+exp(sqrt(2)*(log(time)-beta)/sigma))**2);
			%let log_density = %str(log(sqrt(2)) + sqrt(2)*(log(time)-beta)/sigma - log_sigma - 2*log(max(1+exp(sqrt(2)*(log(time)-beta)/sigma),&log0)));
		%end;
		%if %upcase(&dist.) in WEIBULL GAMMA LNORM LLOGIS %then %do;
	        %let param_anc = sigma;
			%let log_transf_param = sigma;
	    %end;
	%end;

	*generate parameter names after log-transformation;
	%let param_dist = &location &param_anc; *e.g., beta sigma lambda;
	%if &log_transf_param ne %then %do;
		%let param_id = 1;
		%let param_log_dist=; *want: beta log_sigma lambda;
		%let param_log_anc=; *want: log_sigma lambda;
		%do %while (%scan(&param_dist, &param_id) ne );
			%let param_anc_s = %scan(&param_dist, &param_id);
			%let param_label_s = ;
			%if &param_anc_s in &log_transf_param %then %let param_label_s = log_&param_anc_s;
			%else %let param_label_s = &param_anc_s;
			%let param_log_dist = &param_log_dist &param_label_s;
			%if &param_anc_s in &param_anc %then %let param_log_anc = &param_log_anc &param_label_s;
			%let param_id = %eval(&param_id+1);
		%end;
	%end;
	%else %do; *no log-transformation needed;
		%let param_log_dist=&param_dist;
		%let param_log_anc=&param_anc;
	%end;

	%if &dist eq %then %do; *custom distribution;
		*Preparation for likelihood calculation in custom distribution (used in %param_prep);
		*assign values to parameters in the custom model;
		%let param_id = 1;
		%let exec = ;
		%do %while (%scan(&param_dist, &param_id) ne );
			%let param_s = %scan(&param_dist, &param_id);
			%let param_log_s = %scan(&param_log_dist,&param_id);
			%let exec = %str(&exec &param_log_s=param_value[&param_id];);
			%if &log_transf_param ne %then %do;
				%if &param_s in &log_transf_param %then %let exec = %str(&exec &param_s=exp(&param_log_s););
			%end;
			%let param_id = %eval(&param_id+1);
		%end;
	%end;
	
	* check data and delete observations with invalid response;
	data checked;
    	set &data;
    	%if &t2 ne %then %do;
			if (&t1 = . and &t2 = .) or &t1 > &t2 > . or . < &t1 < &log0 or . < &t2 < &log0 then delete;
		%end;
		%else %if &censor ne %then %do;
			if &t1 < &log0 or &censor =. then delete;
		%end;
	run;


	%****** Start IML ******;
    proc iml;


        %****** Calculate log likelihood for each individual *******;

        start loglik_ind(theta) global (times_i, x_i, wt_i, param_vec, cov_index);
            %param_prep();
    		%if %bquote(&survival)= and %bquote(&log_survival)~= %then
        		%let survival = %str(exp(&log_survival));
    		%if %bquote(&density)= and %bquote(&log_density)~= %then
        		%let density = %str(exp(&log_density));
			%if %bquote(&survival)= and %bquote(&log_survival)= and %bquote(&density)~= and %bquote(&hazard)~= %then
        		%let survival = %str((&density) / (&hazard));
    		%if %bquote(&density)= and %bquote(&survival)~= and %bquote(&hazard)~= %then
        		%let density = %str((&survival) * (&hazard));
		    %if %bquote(&hazard)= and %bquote(&density)~= and %bquote(&hazard)~= %then
        		%let hazard = %str((&density) / (&survival));
			%if %bquote(&log_density)= and %bquote(&density)~=%then
        		%let log_density = %str(log(&density));
			%if %bquote(&log_survival)= and %bquote(&survival)~= %then
        		%let log_survival = %str(log(&survival));

            %*right censored time;
            if (times_i[1] ^= .) & (times_i[2] = . ) then do;
                time = times_i[1];
                LL_i = &log_survival*wt_i;
            end;
            %*left censored time;
            if (times_i[1] = .) & (times_i[2] ^= .) then do;
                time = times_i[2];
                LL_i = log(1 - (&survival))*wt_i;
            end;
            %*interval censored time;
            if (times_i[1] ^= .) & (times_i[2] > times_i[1]) then do;
                time = times_i[1];
                temp_val = &survival;
                time = times_i[2];
                y = temp_val - (&survival);
                LL_i = log(y)*wt_i;
            end;
            %*uncensored time;
            if times_i[1] = times_i[2] then do;
                time = times_i[1];
                LL_i = &log_density*wt_i;
            end;
            return(LL_i);
        finish loglik_ind;


    	%****** Calculate overall log likelihood *******;

        start loglik(theta) global (times, x, wt, nobs, times_i, x_i, wt_i);
            LL = 0;
            do i = 1 to nobs;
                times_i = times[i,];
                x_i = x[i,];
                wt_i = wt[i,1];
                LL_i = loglik_ind(theta);
                LL = LL + LL_i;
            end;
            return(LL);
        finish loglik;


		/******* Combine covariates in &covars and &anc and create param_list to store covariates 
		associated with each parameter******/

		param_vec = {}; *parameters in character vector;
		param_vec_log = {}; *parameters (in log scale) in character vector;
		%let param_id = 1;
		%do %while (%scan(&param_dist, &param_id) ne );
			%let param_s = %scan(&param_dist, &param_id);
			%let param_s_log = %scan(&param_log_dist, &param_id);
			param_vec = param_vec || "&param_s";
			param_vec_log = param_vec_log || "&param_s_log";
			%let param_id = %eval(&param_id+1);
		%end;

		covars_com_vec = {};
		%if &covars ne or %bquote(&anc) ne %then %do;
			param_list = ListCreate(param_vec);
			%sublist(param="&location",covars=&covars); *add location parameter to param_list;
			%if %bquote(&anc) ne %then %do;
				%let anc_id = 1;
				%do %while (%qscan(%str(&anc), &anc_id, %str(,)) ne );
					%let anc_s = %qscan(%str(&anc), &anc_id, %str(,));
					%let anc_name = %scan(%str(&anc_s),1,%str(());
					%let anc_cov = %scan(%str(&anc_s),2,%str(());
					%sublist(param="&anc_name",covars=&anc_cov); *add each ancillary parameter to param_list;			
					%let anc_id = %eval(&anc_id+1);
				%end;					
			%end;
		%end;
		
		*get unique values of covars_com_vec but remain the original order;
		covars_combined = unique(covars_com_vec);
		idx = j(ncol(covars_combined),1,0);      
	    do ii = 1 to ncol(covars_combined);
	       idx[ii] = loc(covars_com_vec=covars_combined[ii])[1]; 
	    end;
	    call sort(idx);             
	    covars_combined = T(covars_com_vec[idx]);
		
		covars_combined_string = ""; *string that contains all unique covariates;
		do kk=1 to ncol(covars_combined);
			covar = covars_combined[kk];
			if kk < ncol(covars_combined) then covars_combined_string = covars_combined_string + strip(covar) + " ";
			else covars_combined_string = covars_combined_string + strip(covar);
		end;


		%******* Stratification set-up ******;

		%if &strata ne %then %do; *with stratification;
            use checked;
            	read all var {&strata} into factors;
            	read all var _ALL_;
            close checked;
            nfactor = ncol(factors);
    		if type(factors)='N' then factor = strip(putn(factors[,1],"BEST6."));
			else factor = factors[,1];
            if nfactor > 1 then
                do fac=2 to nfactor;
                    factor = catx(", ", factor, factors[,fac]);
                end;
			fac_var = {&t1 &t2 &censor} || covars_combined || {&weight factor};
            create data_fac var fac_var; *create a new dataset called data_fac to store strata information;
            append;
            levels = unique(factor);
            call symputx("nlevels", ncol(levels));
            %let wh_statement = %str(where(factor=st));
            %let dataset = data_fac;
        %end;
        %else %do; *without stratification;
            levels=.;
            %let nlevels=1;
            %let wh_statement =;
            %let dataset = checked;
        %end;
		nobs_level = j(%eval(&nlevels), 1, 0);
		est_level = {};
		sd_level = {};


		%******* Fitting model within each stratum ******;
		
        %do k=1 %to &nlevels;
            st = levels[&k];
			call symputx("st", st);  
			%if &strata ne %then %let wh_st_x = %nrstr(where=(factor='&st'));
			%else %let wh_st_x =;

			%******* Load response and weights ******;

            use &dataset &wh_statement;
				read all var {&t1} into t1;
				nobs = nrow(t1);
				nobs_level[&k] = nobs;
				%if &t2 eq and &censor ne %then %do; *if using t1/censor, convert it to t1/t2;
					times = t1 || j(nobs,1,0);
					read all var {&censor} into censor;
					do i = 1 to nobs;
		        		if censor[i] = %eval(&censval) then do; times[i,2] = .; end;
		        		else do; times[i,2] = times[i,1]; end;
					end;
		    	%end;
				%else %do; read all var {&t1 &t2} into times; %end;
                if %eval(&weight=) then wt = j(nobs,1,1);
                else read all var {&weight} into wt;
            close &dataset;
			
			%******* Load and process covariates ******;

			col_one = j(nobs,1,1);
			cov_index = ListCreate(param_vec); *list that stores column indices for each parameter;
            %if &covars eq and %bquote(&anc) eq %then %do; *no covariates in the model;
                x=col_one; *column of one;
                if %eval(&param_anc=) then param_log={"&location"}`;
                else param_log = {&param_log_dist}`;
                if %eval(&param_anc=) then param={"&location"}`;
                else param = {&param_dist}`;
				do param_id=1 to ncol(param_vec);
					call ListSetItem(cov_index,param_id,{1}); *index={1} only intercept;
				end;
            %end; 
            %else %do; *there are covariates in the model;
                %class_cov_prep(dataset=&dataset,where_st=%str(&wh_st_x)); *prepare continuous data matrix (including dummies);
                x = col_one||x_num; *add column of one;
				%if %bquote(&anc) eq %then %do;     *no ancillary regression;
					call ListSetItem(cov_index, "&location", do(1,ncov_com+1,1)); *use all covariates;
					cov_name = {};
					do cov_i=1 to ncov_com;
						covar = ListGetItem(param_list,"&location")[cov_i];	
						cov_name = cov_name // ListGetItem(ListGetItem(cov_name_dic,covar),"covarLabel");
					end;
					*create labels for output;
	                if %eval(&param_anc=) then param_log={"Intercept"} // cov_name;
	                else param_log = {"Intercept"} // cov_name //{&param_log_anc}`;
	                if %eval(&param_anc=) then param="Intercept" // cov_name;
	                else param = {"Intercept"}// cov_name //{&param_anc}`;
					do param_id=2 to ncol(param_vec);
						param_s = param_vec[param_id];
						call ListSetItem(cov_index, param_s, {1}); *index={1} for other parameters;
					end;
				%end;
				%else %do;     *ancillary regression;
					param = {};
					param_log = {};
					do param_id=1 to ncol(param_vec); *for each parameter;
						param_s = param_vec[param_id]; *parameter name;
						param_s_log = param_vec_log[param_id];
						covar_vec = ListGetItem(param_list,param_s);
						if IsEmpty(covar_vec) then do; *no covariate for this parameter;
							param = param // param_s;
							param_log = param_log // param_s_log;
							index = {1};
						end;
						else do cov_i=1 to ncol(covar_vec); *there are covariates for this parameter, for each covariate;
							covar = covar_vec[cov_i]; *covariate name;
							if cov_i=1 then do; *intercept;
								param = param // catx(": ", param_s,"intercept");
								param_log = param_log // catx(": ", param_s_log,"intercept");
								index = {1};
							end;
							cov_list = ListGetItem(cov_name_dic,covar); *sublist for this covariate;
							param = param // catx(": ", param_s, ListGetItem(cov_list,"covarLabel"));
							param_log = param_log // catx(": ", param_s_log, ListGetItem(cov_list,"covarLabel"));
							index = index || (ListGetItem(cov_list,"covarIndex")+1); * 1 is for intercept;
						end;
						call ListSetItem(cov_index, param_s, index);
					end;
				%end;
            %end;		
			nparams = nrow(param);


			%******* Initial values ******;

			%if %bquote(&init) eq %then %do; * if no initial values supplied;
				theta0 = j(nparams, 1, 0); 
				t = j(nobs,1,0);
				log_t = j(nobs,1,0);
				do i = 1 to nobs;
					if (times[i,1] ^= .) & (times[i,2] > times[i,1]) then t[i]=(times[i,1]+times[i,2])/2;
					else if (times[1] = .) & (times[2] ^= .) then t[i] = times[i,2];
					else t[i] = times[i,1];
					if t[i]>0 then log_t[i] = log(t[i]);
				end;
				theta0[1] = log_t[:];	*initial value for beta;
				%if &dist ne %then %do;     *built-in distribution;				
					%if %upcase(&dist.) in GENGAMMA GENF GENF_ORIG WEIBULL GAMMA LNORM LLOGIS %then %do;
						ncovs_beta = ncol(ListGetItem(cov_index, "beta"));
						ncovs_sigma = ncol(ListGetItem(cov_index, "sigma"));
						if ncovs_sigma>1 then theta0[ncovs_beta+1] = log(std(log_t));     *initial value for sigma;
						else theta0[ncovs_beta+1] = std(log_t);
						%if %upcase(&dist.) in GENGAMMA GENF GENF_ORIG %then %do;
							theta0[ncovs_beta+ncovs_sigma+1] = 0.5;     *initial value for lambda / q;
						%end;
						%if %upcase(&dist.) = GENF %then %do;
							ncovs_q = ncol(ListGetItem(cov_index, "q"));
							ncovs_p = ncol(ListGetItem(cov_index, "p"));
							if ncovs_sigma=1 then theta0[ncovs_beta+ncovs_sigma+ncovs_q+1] = 1;     *initial value for p;
						%end;
						%if %upcase(&dist.) = GENF_ORIG %then %do;
							ncovs_m1 = ncol(ListGetItem(cov_index, "m1"));
							if ncovs_m1=1 then theta0[ncovs_beta+ncovs_sigma+1] = 1;     *initial value for m1;
							ncovs_m2 = ncol(ListGetItem(cov_index, "m2"));
							if ncovs_m2=1 then theta0[ncovs_beta+ncovs_sigma+ncovs_m1+1] = 1;     *initial value for m2;
						%end;
					%end;
					%if %upcase(&dist.)=GOMPERTZ %then %do;
						ncovs_beta = ncol(ListGetItem(cov_index, "beta"));
						theta0[ncovs_beta+1] = 0.5;
					%end;
				%end;
				%else %if &log_transf_param ne %then %do;     *custom distribution;
					%let log_transf_id = 1;
					%do %while (%scan(&log_transf_param, &log_transf_id) ne ); *for covariates that need log-transformation;
						%let param_s = %scan(&log_transf_param, &log_transf_id);
						param_s_string = "&param_s";
						if ncol(ListGetItem(cov_index, param_s_string))=1 then do;
							index = loc(upcase(param)=upcase(param_s_string));
							theta0[index] = 1; *set initial value of 1;
						end;
						%let log_transf_id = %eval(&log_transf_id+1);
					%end;
				%end;
			%end;
            %else %do; theta0 = &init.`; %end;     *user-provided initial values;
			theta0_transf = theta0;     *initial values in the original scale;

			* transform from original scale to log scale;
			%let log_transf_id = 1;
			%do %while (%scan(&log_transf_param, &log_transf_id) ne );
				%let param_s = %scan(&log_transf_param, &log_transf_id);
				param_s_string = "&param_s";
				if ncol(ListGetItem(cov_index, param_s_string))=1 then do;
					index = loc(upcase(param)=upcase(param_s_string));
					theta0[index] = log(theta0[index]);
				end;
				%let log_transf_id = %eval(&log_transf_id+1);
			%end;


			%****** Call optimization function ******;

            a = &lower;
            b = &upper;
            con = a//b;
            con = con[, 1:nparams];
            optn = {1,%eval(&nlp_print)};
            call &optim_method(rc,thetares,"loglik",theta0,optn,con);
            thetaopt_log = thetares`; *optimal theta in log scale for positive parameters;
            maxll = loglik(thetaopt_log); *maximum log-likelihood;
			%if &dist ne %then %do;
				%if %upcase(&dist.) in GENGAMMA GENF GENF_ORIG WEIBULL GAMMA LLOGIS %then %do;
					do i=1 to nobs; *adjust log-likelihood;
						if times[i,1] = times[i,2] then maxll = maxll - log(times[i,1])*wt[i];
					end;
	            %end;
			%end;
			aic = -2*maxll+2*nparams;
			bic = -2*maxll+log(nobs)*nparams;
			* transform from log scale to original scale;
			thetaopt = thetaopt_log;
			%let log_transf_id = 1;
			%do %while (%scan(&log_transf_param, &log_transf_id) ne );
				%let param_s = %scan(&log_transf_param, &log_transf_id);
				param_s_string = "&param_s";
				if ncol(ListGetItem(cov_index, param_s_string))=1 then do;
					index = loc(upcase(param)=upcase(param_s_string));
					thetaopt[index] = exp(thetaopt_log[index]);
				end;
				%let log_transf_id = %eval(&log_transf_id+1);
			%end;
			est_level = est_level//thetaopt`;


            %****** Inference ******;

            call nlpfdd(LL,deriv,h,"loglik", thetaopt_log);
            cov = inv(-h);
            %let label_se = S.E.;
            %if %upcase(&robust)=YES %then %do;     *Sandwich estimator;
                B = j(nparams, nparams, 0);
                do i = 1 to nobs;
                    times_i = times[i,];
                    x_i = x[i,];
                    call nlpfdd(LL,deriv,h,"loglik_ind", thetares);
                    U_i = deriv`*deriv;
                    B = B+U_i; 
                end;
                cov = cov * B * cov;
                %let label_se = %str(Robust S.E.);;
            %end;
            setheta_log = sqrt(vecdiag(cov)<>0);

			* transform from log scale to original scale;
			setheta = setheta_log;
			%let log_transf_id = 1;
			%do %while (%scan(&log_transf_param, &log_transf_id) ne );
				%let param_s = %scan(&log_transf_param, &log_transf_id);
				param_s_string = "&param_s";
				if ncol(ListGetItem(cov_index, param_s_string))=1 then do;
					index = loc(upcase(param)=upcase(param_s_string));
					setheta[index] = setheta_log[index]*thetaopt[index];
				end;
				%let log_transf_id = %eval(&log_transf_id+1);
			%end;

			* confidence intervals;
            ci_l_log = thetaopt_log - quantile('normal',1-%sysevalf(&alpha/2)) * setheta_log;
            ci_u_log = thetaopt_log + quantile('normal',1-%sysevalf(&alpha/2)) * setheta_log;
			ci_l = ci_l_log;
            ci_u = ci_u_log;
			%let log_transf_id = 1;
			%do %while (%scan(&log_transf_param, &log_transf_id) ne );
				%let param_s = %scan(&log_transf_param, &log_transf_id);
				param_s_string = "&param_s";
				if ncol(ListGetItem(cov_index, param_s_string))=1 then do;
					index = loc(upcase(param)=upcase(param_s_string));
					ci_l[index] = exp(ci_l_log[index]);
					ci_u[index] = exp(ci_u_log[index]);
				end;
				%let log_transf_id = %eval(&log_transf_id+1);
			%end;
			sd_level = sd_level//setheta`;

            t_value_log = thetaopt_log / setheta_log;
            p_value_log = 2*(1-cdf('normal',abs(t_value_log)));
			t_value = thetaopt / setheta;
            p_value = 2*(1-cdf('normal',abs(t_value)));


            %****** prediction ******;

			%if &pred eq %then %goto output;
			%if &covars ne or %bquote(&anc) ne %then %do; *if there are covariates, prepare the covariates to continuous matrix;
				%class_cov_prep(dataset=&pred);
			%end;
            use &pred;
                read all var {time} into pred_time;
            close &pred;
            n_pred = nrow(pred_time);
            survival_pred = j(n_pred, 1, 0);
            se_surv_pred = j(n_pred, 1, 0);
            hazard_pred = j(n_pred, 1, 0);
            se_haz_pred = j(n_pred, 1, 0);
			col_one_pred = j(n_pred, 1, 1);
            %if &covars eq and %bquote(&anc) eq %then %do; x_pred = col_one_pred; %end;
            %else %do; x_pred = col_one_pred||x_num; %end;

			start surv(theta) global (time, x_i, param_vec, cov_index); *predicted survival;
				%param_prep();
                surv_pred = &survival;
                return(surv_pred);
            finish surv;

			start haz(theta) global (time, x_i, param_vec, cov_index); *predicted hazard;
				%param_prep();
				%let hazard = %str((&density) / (&survival));
                hazard_pred = &hazard;
                return(hazard_pred);
            finish haz;

            do i = 1 to n_pred;
				time = pred_time[i];
                x_i = x_pred[i,];
                
                call nlpfdd(s, deriv_s, hess_s, "surv",thetaopt_log);
                survival_pred[i] = s;
                se_surv_pred[i] = sqrt(deriv_s*cov*deriv_s`<>0);

                call nlpfdd(hz, deriv_hz, hess_hz, "haz",thetaopt_log);
                hazard_pred[i] = hz;
                se_haz_pred[i] = sqrt(deriv_hz*cov*deriv_hz`<>0);
            end;
			pred_r = TableCreateFromDataSet("Work", "&pred", "keep = time "+covars_combined_string);
			call TableAddVar(pred_r, "Survival", survival_pred);
			call TableAddVar(pred_r, "Survival SE", se_surv_pred);
			call TableAddVar(pred_r, "Hazard", hazard_pred);
			call TableAddVar(pred_r, "Hazard SE", se_haz_pred);
			call TableSetVarFormat(pred_r, {"Survival" "Survival SE" "Hazard" "Hazard SE"}, {"7.3" "7.3" "7.3" "7.3"});


			%******Predicted survival and hazard curves ******;

			%if &pred_max_time= %then %do; t = do(max(times)/100, max(times),max(times)/100); %end;
			%else %do; t = do(%sysevalf(&pred_max_time/100), %eval(&pred_max_time),%sysevalf(&pred_max_time/100)); %end;
			
			* remove duplicated rows in x_pred;
			dupRows = j(n_pred, 1, 0);
			do ii = 1 to n_pred-1;
      			if dupRows[ii] = 0 then do; 
         			r = ii+1:n_pred;           
         			M = x_pred[r,]-x_pred[ii,];   
         			b = M[,##];       
        			if ncol(loc(b=0)) > 0 then dupRows[r[loc(b=0)]] = 1;
      			end;      
   			end;
			uniqRows = loc(dupRows=0);
			x_pred = x_pred[uniqRows, ];
			n_pred = nrow(x_pred);

			%if &covars eq and %bquote(&anc) eq %then %do;	ncov_level = 1; %end;
			%else %do; ncov_level = n_pred; %end;
			if %eval(&k=1) then do;
				ss_res = j(ncov_level*ncol(t),6,1);     *time, predicted survival, lower CI, upper CI, covariate_id, stratum_id;
				hh_res = j(ncov_level*ncol(t),6,1);     *time, predicted hazard, lower CI, upper CI, covariate_id, stratum_id;
				ss_res[,1] = repeat(t`,ncov_level,1);
				hh_res[,1] = repeat(t`,ncov_level,1);
				ss_res[,6] = %eval(&k);
				hh_res[,6] = %eval(&k);
				group = {};
				group_id = {};
			end;
			else do;
				ap = repeat(t`,ncov_level,1) || j(ncov_level*ncol(t),4,1) || j(ncov_level*ncol(t),1,%eval(&k));     *append;
				ss_res = ss_res // ap; 
				hh_res = hh_res // ap;			
			end;	
			%if &covars eq and %bquote(&anc) eq %then %do; *no covariates;
				x_i = 1;
				do j = 1 to ncol(t);
					time = t[j];			
					qu = quantile('normal',1-%sysevalf(&alpha/2));
					call nlpfdd(s, deriv_s, hess_s, "surv",thetaopt_log);
	                st_se = sqrt(deriv_s*cov*deriv_s`<>0);
	                ss_res[nrow(ss_res)-ncol(t)+j,2:4] = s||(s - qu* st_se)||(s + qu* st_se);
	                call nlpfdd(hz, deriv_hz, hess_hz, "haz",thetaopt_log);
	                ht_se = sqrt(deriv_hz*cov*deriv_hz`<>0);		
	                hh_res[nrow(hh_res)-ncol(t)+j,2:4] = hz||(hz - qu* ht_se)||(hz + qu* ht_se);
				end;
				%if &strata eq %then %do; %let grp = ; %end;
				%else %do;
					group = group // st;
					group_id = group_id // catx(',', 1, &k);
					%let grp = %str(group=group);
				%end;
			%end;
			%else %do; *with covariates;
				do i = 1 to n_pred;
					ss_res[(nrow(ss_res)-(n_pred-i+1)*ncol(t)+1):(nrow(ss_res)-(n_pred-i)*ncol(t)),5]=i;
					hh_res[(nrow(ss_res)-(n_pred-i+1)*ncol(t)+1):(nrow(ss_res)-(n_pred-i)*ncol(t)),5]=i;
	               	x_i = x_pred[i,];
					*create labels for groups;	
					string = "";
					do kk=1 to ncol(covars_combined);
						covar = covars_combined[kk];
						xc = TableGetVarData(x_table, covar);
						xc = xc[i];
						if type(xc)='N' then xc = putn(xc,"BEST6.");
						if kk < ncol(covars_combined) then string = string + covar + "=" + strip(xc) + ", ";
						else string = string + covar + "=" + strip(xc);
					end;
					if %eval(&strata eq) then group = group // string;
					else do;
						if type(st)='N' then st = putn(st,"BEST6.");
						group = group // (string + ", " + "&strata=" + st);
					end;
					group_id = group_id // catx(',', i, &k);
					do j = 1 to ncol(t);
						time = t[j];			
						qu = quantile('normal',1-%sysevalf(&alpha/2));
						call nlpfdd(s, deriv_s, hess_s, "surv",thetaopt_log);
		                st_se = sqrt(deriv_s*cov*deriv_s`<>0);		
		                ss_res[nrow(ss_res)-(n_pred-i+1)*ncol(t)+j,2:4] = s||(s - qu* st_se)||(s + qu* st_se);
		                call nlpfdd(hz, deriv_hz, hess_hz, "haz",thetaopt_log);
		                ht_se = sqrt(deriv_hz*cov*deriv_hz`<>0);		
		                hh_res[nrow(ss_res)-(n_pred-i+1)*ncol(t)+j,2:4] = hz||(hz - qu* ht_se)||(hz + qu* ht_se);
					end;
				end;
				%let grp = %str(group=group);
			%end;


            %****** output ******;

			%output: 			
			%if %upcase(&res_print)=NO %then %goto fin;
			%if &k=1 %then %do;
				*Specifications;
				tb_spec = TableCreate('Feature', {"Data Set","Distribution","Optimization Method"});
				%if &dist ne %then %do;	dist_output ="&dist"; %end;
				%else %do; dist_output = "Custom"; %end;
				specs = {"&data"}//dist_output//{"&optim_method"};
				call TableAddVar(tb_spec, {"Value"},specs);
				call TablePrint(tb_spec) ID='Feature' label='Specifications';
				*Initial values;
				tb_init = TableCreate('Variable', param);
				call TableAddVar(tb_init, {"Value"},theta0_transf);
				call TableSetVarFormat(tb_init, {"Value"}, {'D11.3'});
				call TablePrint(tb_init) ID='Variable' label='Initial Parameters';
			%end;
			*Strata info;
            if %eval(&strata ne ) then do;
               strata_info = cat('Stratum ',&k,': ',"&strata"," = ",st, ", n = ",nobs_level[&k]);
               print strata_info[label=""];
            end;
			*Estimation;
            if rc < 0 then print "NOTE: Iteration failed to converge.  Estimates are unreliable.";
			tb_est = TableCreate('Variable', param);
			call TableAddVar(tb_est, {"Estimate" "&label_se" "95% CI Lower" "95% CI Upper" "t Value" "Pr>|t|"},
                        thetaopt || setheta || ci_l || ci_u || t_value || p_value);
			call TableSetVarFormat(tb_est, {"Estimate" "&label_se" "95% CI Lower" "95% CI Upper" "t Value" "Pr>|t|"},
                            {'D11.3' 'D11.3' 'D11.3' 'D11.3' '7.2' 'PVALUE6.4'});
			call TablePrint(tb_est) ID='Variable' label='Parameter Estimates';
			if %eval(%upcase(&log_result)=YES) then do; *results in the log scale;
				tb_est_log = TableCreate('Variable', param_log);
				call TableAddVar(tb_est_log, {"Estimate" "&label_se" "95% CI Lower" "95% CI Upper" "t Value" "Pr>|t|"},
	                        thetaopt_log || setheta_log || ci_l_log || ci_u_log || t_value_log || p_value_log);
				call TableSetVarFormat(tb_est_log, {"Estimate" "&label_se" "95% CI Lower" "95% CI Upper" "t Value" "Pr>|t|"},
	                            {'D11.3' 'D11.3' 'D11.3' 'D11.3' '7.2' 'PVALUE6.4'});
				call TablePrint(tb_est_log) ID='Variable' label='Parameter Estimates';
			end;
			*Log-likelihood;
			print maxll[l="Log-likelihood" F=7.3] aic[l="AIC" F=7.3] bic[l="BIC" F=7.3];
			*Prediction;
			if %eval(&pred ne ) then run TablePrint(pred_r) label="Prediction";
        %end;


		%****** Combined results from various strata ******;

		%if &strata ne %then %do;
			wt_level = nobs_level/sum(nobs_level); *weight;
			thetaopt_level = est_level`*wt_level;
			setheta_level = sqrt((sd_level#sd_level)`*(wt_level#wt_level)<>0);
			t_value_level = thetaopt_level / setheta_level;
            p_value_level = 2*(1-cdf('normal',abs(t_value_level)));
            ci_l_level = thetaopt_level - quantile('normal',1-%sysevalf(&alpha/2)) * setheta_level;
            ci_u_level = thetaopt_level + quantile('normal',1-%sysevalf(&alpha/2)) * setheta_level;
			tb_stata = TableCreate('Variable', param);
			call TableAddVar(tb_stata, {"Estimate" "&label_se" "95% CI Lower" "95% CI Upper" "t Value" "Pr>|t|"},
                        thetaopt_level || setheta_level || ci_l_level || ci_u_level || t_value_level || p_value_level);
			call TableSetVarFormat(tb_stata, {"Estimate" "&label_se" "95% CI Lower" "95% CI Upper" "t Value" "Pr>|t|"},
                            {'D11.3' 'D11.3' 'D11.3' 'D11.3' '7.2' 'PVALUE6.4'});
			call TablePrint(tb_stata) ID='Variable' label='Combined results using sample size weights';
		%end;


		%****** Create prediction dataset for figures ******;

		%if &pred ne %then %do;
			create pred_surv from ss_res[colname = {'time' 'survival' 'lower_ci' 'upper_ci' 'covariate' 'Strata'}];
	        append from ss_res;
			create pred_haz from hh_res[colname = {'time' 'hazard' 'lower_ci' 'upper_ci' 'covariate' 'Strata'}];
           	append from hh_res;
			%if &covars ne or %bquote(&anc) ne or &strata ne %then %do;
				create cov_name_pred var {"group_id" "group"};
				append;
			%end;
		%end;
    quit;


	%****** Predicted survival/hazard figures ******;

	%if &pred ne %then %do;
		%if &covars ne or %bquote(&anc) ne or &strata ne %then %do;
			data pred_surv;
				set pred_surv;
				group_id = catx(',', covariate, Strata); 
			run;
			data pred_haz;
				set pred_haz;
				group_id = catx(',', covariate, Strata);
			run;
			proc sort data = cov_name_pred; by group_id; run;
			proc sort data = pred_surv; by group_id; run;
			proc sort data = pred_haz; by group_id; run;
			data pred_surv;	merge pred_surv cov_name_pred; by group_id;	run;
			data pred_haz;	merge pred_haz cov_name_pred; by group_id;	run;
		%end;
		proc sgplot data=pred_surv;
   			series x=time y=survival / &grp;
			%if %upcase(&pred_plot_cl)=YES %then %do; 
				band x=time lower=lower_ci upper=upper_ci /&grp transparency=0.5;
			%end;
   			yaxis label="Predicted Survival";
		run;
		proc sgplot data=pred_haz;
   			series x=time y=hazard / &grp;
			%if %upcase(&pred_plot_cl)=YES %then %do; 
				band x=time lower=lower_ci upper=upper_ci /&grp transparency=0.5;
			%end;
   			yaxis label="Predicted Hazard";
		run;
	%end;

	%fin:
%mend paramsurv;






