*########################################
*######################################## Program Simulations
*########################################

capture program drop AR_sim 
capture program drop MA_sim 
capture program drop ARMA_sim 

*##################################################
*################################################## AR simulation
*##################################################

program AR_sim
syntax  ,[LENgth_vector(numlist int max=1) LAG_p(numlist int max=1) Mean_ra(numlist  max=1) SD_ra(numlist  max=1)  COEFficient(numlist max=1) INt_value(numlist  max=1) ]
if "`Int_Value'"==""{
display "The order of arguments are length_vector lag_p mean_ra sd_ra coefficient Initial_Value"
}

if "`length_vector'"==""{
local length_vector = _N
}

if "`lag_p'"==""{
local lag_p = 1
}

if "`mean_ra'"==""{
local mean_ra = float(0)
}

if "`sd_ra'"==""{
local sd_ra = float(1)
}

if "`coefficient'"==""{
local coefficient = float(0.5)
}

if "`int_value'"==""{
local int_value = 0
}

** Parameters Input:
local largo_Vect = `length_vector'
local rezagos_p = `lag_p'
local media = `mean_ra'
local des_E = `sd_ra'
local coeficiente = `coefficient'
local Int_Value = `int_value'

*##################################################
*################################################## Auto regressive simulation
*##################################################

** Empty Matrix
mat M = J( `largo_Vect' ,1,0 )  // Vector for save output values
mat E_V = J( `largo_Vect' ,1,0 ) // Vector for errors

** Forcast first Part
if `rezagos_p' == 1{
						mat M[1,1] = `Int_Value' //Put initial value in Matrix
						}
						else{
						mat M[1,1] = `Int_Value' //Put initial value in Matrix
						forvalues i=2/`rezagos_p'{
							local i_1 = `i' - 1						
							mat A = J(1,`i_1', `coeficiente')  		
							forvalues j = 1/`i_1'{
									local expon = `i' - `j' + 1
									mat A[1, `j'] = A[1, `j']^`expon'
							}
							local i_1 = `i' - 1	
							if `i_1'== 1{
							mat U = M[1,1]
							}
							else{
							mat U = M[1..`i_1',1]
							}
							local error_t = rnormal(`media', `des_E')
							mat uu_t = A*U
							local u_t = uu_t[1,1] + `error_t'
							mat M[`i',1] = `u_t'
						}
}
** Forcast second Part
mat A = J(1, `rezagos_p', `coeficiente' )
forvalues j = 1/`rezagos_p'{
					local expon = `rezagos_p' - `j' + 1
					mat A[1, `j'] = A[1, `j']^`expon'
					}

local i_b = `rezagos_p' + 1
forvalues i= `i_b'/`largo_Vect'{
local error_t = rnormal(`media', `des_E')
local i_1 = `i' - 1
local i_p = `i' - `rezagos_p' 
mat U = M[`i_p'..`i_1', 1]
mat uu_t = A*U
local u_t = uu_t[1,1] + `error_t'
mat M[`i',1]=`u_t'
}
svmat  M	
rename M1 AR_`lag_p'
end

*##################################################
*################################################## MA simulation
*##################################################

program MA_sim
syntax  ,[LENgth_vector(numlist int max=1) LAG_q(numlist int max=1) Mean_ma(numlist  max=1) SD_ma(numlist  max=1)  COEFficient(numlist max=1) INt_value(numlist  max=1) ]
if "`Int_Value'"==""{
display "The order of arguments are length_vector lag_q mean_ma sd_ma coefficient Initial_Value"
}

if "`length_vector'"==""{
local length_vector = _N
}

if "`lag_q'"==""{
local lag_q = 1
}

if "`mean_ma'"==""{
local mean_ma = 0
}

if "`sd_ma'"==""{
local sd_ma = 1
}

if "`coefficient'"==""{
local coefficient = 0.5
}

if "`Int_Value'"==""{
local Int_Value = 1
}
** Parameters Input:
local largo_Vect = `length_vector'
local rezagos_p = `lag_q'
local media = `mean_ma'
local des_E = `sd_ma'
local coeficiente = `coefficient'
local Int_Value = `Int_Value'

*##################################################
*################################################## MA
*##################################################

** Empty Matrix
mat M = J( `largo_Vect' ,1,0 )  // Vector for save output values
mat E_V = J( `largo_Vect' ,1,0 ) // Vector for errors

** Forcast first Part
if `rezagos_p' == 1{
						mat M[1,1] = `Int_Value' //Put initial value in Matrix
						}
						else{
						mat M[1,1] = `Int_Value' //Put initial value in Matrix
						forvalues i=2/`rezagos_p'{
							local i_1 = `i' - 1						
							mat A = J(1,`i_1', `coeficiente')  		
							forvalues j = 1/`i_1'{
									local expon = `i' - `j' + 1
									mat A[1, `j'] = A[1, `j']^`expon'
							}
							local i_1 = `i' - 1	
							if `i_1'== 1{
							mat U = E_V[1,1]
							}
							else{
							mat U = E_V[1..`i_1',1]
							}
							mat E_V[`i',1] = rnormal(`media', `des_E')
							mat uu_t = A*U
							local u_t = uu_t[1,1] + E_V[`i',1]
							mat M[`i',1] = `u_t'
						}
}
** Forcast second Part
mat A = J(1, `rezagos_p', `coeficiente' )
forvalues j = 1/`rezagos_p'{
					local expon = `rezagos_p' - `j' + 1
					mat A[1, `j'] = A[1, `j']^`expon'
					}

local i_b = `rezagos_p' + 1
forvalues i = `i_b'/`largo_Vect'{
mat E_V[`i',1] = rnormal(`media', `des_E')
local i_1 = `i' - 1
local i_p = `i' - `rezagos_p' 
mat U = E_V[`i_p'..`i_1', 1]
mat uu_t = A*U
local u_t = uu_t[1,1] + E_V[`i',1]
mat M[`i',1]=`u_t'
}
svmat  M	
rename M1 MA_`lag_q'
end



*##################################################
*################################################## ARMA simulation
*##################################################
program ARMA_sim
syntax  ,[LENgth_vector(numlist int max=1) LAG_P(numlist int max=1) LAG_Q(numlist int max=1) Mean_arma(numlist  max=1) SD_arma(numlist  max=1)  COEFFICIENT_AR(numlist max=1) COEFFICIENT_MA(numlist max=1) INt_value(numlist  max=1) ]
args length_vector lag_p lag_q mean_arma sd_arma coefficient_ar coefficient_ma Int_Value
if "`Int_Value'"==""{
display "The order of arguments are length_vector lag_p lag_q mean_arma sd_arma coefficient_ar coefficient_ma  Initial_Value"
}
if "`length_vector'"==""{
local length_vector = _N
}
if "`lag_p'"==""{
local lag_p = 1
}
if "`lag_q'"==""{
local lag_q = 1
}
if "`mean_arma'"==""{
local mean_arma = 0
}
if "`sd_arma'"==""{
local sd_arma = 1
}
if "`coefficient_ar'"==""{
local coefficient_ar = 0.5
}
if "`coefficient_ma'"==""{
local coefficient_ma = 0.2
}
if "`Int_Value'"==""{
local Int_Value = 1
}
** Parameters Input:
local largo_Vect = `length_vector'
local rezagos_p = `lag_p'
local rezagos_q = `lag_q'
local media_ma = `mean_arma'
local des_E_ma = `sd_arma'
local coeficiente_p = `coefficient_ar'
local coeficiente_q = `coefficient_ma'
local Int_Value_ra = `Int_Value'



*##################################################
*################################################## MA
*##################################################

** Empty Matrix
mat M = J( `largo_Vect' ,1,0 )  // Vector for save output values
mat E_V = J( `largo_Vect' ,1,0 ) // Vector for errors

if `rezagos_p' == 1{
						mat M[1,1] = `Int_Value_ra' //Put initial value in Matrix
						}
						else{
						mat M[1,1] = `Int_Value_ra' 						//Put initial value in Matrix
						forvalues i = 2/`rezagos_p'{   					// Begin loop for lag > 2 
							local i_1 = `i' - 1						
							mat A = J(1,`i_1', `coeficiente_p')  		// Empty Vector
							forvalues j = 1/`i_1'{
									local expon = `i' - `j' + 1			// lag coefficient Exponent
									mat A[1, `j'] = A[1, `j']^`expon'
							}
						
						// Section for MA variation
						if `i' <= `rezagos_q'{
							local i_1 = `i' - 1
							mat B = J(1,`i_1', `coeficiente_q' )
							forvalues j = 1/`i_1'{
								local exponq = `i' - `j' + 1 
								mat B[1,`j'] = B[1, `j']^`exponq'
							}
							mat U = E_V[1..`i_1',1]
						}
						else{
							mat B=J(1, `rezagos_q', `coeficiente_q')
							forvalues j = 1/`rezagos_q'{
								local exponq = `rezagos_q' - `j' + 1 
								mat B[1,`j'] = B[1, `j']^`exponq'
							}
							mat U = E_V[1..`rezagos_q',1]
						}
						mat E_V[`i',1] = rnormal( `media_ma', `des_E_ma'  )
						mat ee_t = B*U
						local error_t = ee_t[1,1] + E_V[`i',1]						
						local i_1 = `i' - 1	
							if `i_1'== 1{
							mat ALPHA = M[1,1]
							}
							else{
							mat ALPHA = M[1..`i_1',1]
							}
						mat uu_t = A*ALPHA
						local u_t = uu_t[1,1] + `error_t'
						mat M[`i',1] = `u_t'
						}
						}

mat A = J(1,`rezagos_p', `coeficiente_p')  		// Empty Vector
forvalues j = 1/`rezagos_p'{
	local expon = `rezagos_p' - `j' + 1			// lag coefficient Exponent
	mat A[1, `j'] = A[1, `j']^`expon'
}						
				
						
local i_b = `rezagos_p' + 1	
 
forvalues i = `i_b'/`largo_Vect'{
if  `i' <= `rezagos_q' {
local i_1 = `i'-1
mat B=J(1, `i_1', `coeficiente_q')
forvalues j = 1/`i_1'{
	local exponq = `i' - `j'
	mat B[1,`j'] = B[1, `j']^`exponq'
}
mat U = E_V[1..`i_1',1]
}		
else{
mat B=J(1, `rezagos_q', `coeficiente_q')
forvalues j = 1/`rezagos_q'{
	local exponq = `rezagos_q' - `j' + 1 
	mat B[1,`j'] = B[1, `j']^`exponq'
}
mat U = E_V[1..`rezagos_q',1]
}

mat E_V[`i',1]=rnormal( `media_ma',`des_E_ma' )
mat ee_t = B*U
local error_t = ee_t[1,1] + E_V[`i',1]
local i_1 = `i' - 1
local i_p = `i' - `rezagos_p' 
mat ALPHA = M[`i_p'..`i_1',1]
mat uu_t = A*ALPHA
local u_t = uu_t[1,1] + `error_t'
mat M[`i',1] = `u_t' 
}

svmat  M	
rename M1 ARMA_`lag_p'_`lag_q'

end













