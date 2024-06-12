
;; 1. Based on: run_onco_001
;; 2. Description: Modelling combination effect of Erlotinib and Dactolisib with ATUX
;; x1. Author: My-Luong.Vuong
;; 3. Label: CELLINE_DRUG_ADDITIVE_ADDITIVECONC as ID
;; 4. Author: My-Luong Vuong

$SIZES LTH=9999 LVR=9999
       DIMCNS=9999 DIMTMP=9999 DIMNEW=9999
;        DIMVRB=9999 LIM1=9999 MAXFCN=9999 NO=9999 MAXPTHETA=9999
$PROBLEM Quantitative drug sensitivity scoring of 2 coumpounds in 2 strains

$INPUT CELLINE         ; 2 cell lines: WT, R183W
       DRUG            ; Dactolisib, Erlotinib
       ADDITIVE        ; ATUX
       ADDITIVECONC    ; 0, 3, 4.5, 6 µM
       DAY             ; inter-day repeats
       REP             ; within-day repeats
       DRUGCONC        ;  0.00000, 0.00076, 0.00229, 0.00686, 0.02058, 0.06173, 0.18519, 0.55556, 1.66667, 0.02667, 0.08002, 0.24005, 0.72016, 2.16049, 6.48148, 19.44444, 58.33333 
       DV           
       DROP ; ID option 1 --- CELLINE_DRUG_ADDITIVE
       ID   ; ID option 2 --- CELLINE_DRUG_ADDITIVE_ADDITIVECONC ; --- Estimate parameters for each value of CELLINE_DRUG_ADDITIVE_ADDITIVECONC = combination of the covariates, with REPEAT in the IOV/RUV and REP in the RUV, or ... define CELLINE_DRUG_ADDITIVE_ADDITIVECONC as ID ...
       DROP ; ID option 3 --- CELLINE_DRUG_ADDITIVE_ADDITIVECONC_DAY ; CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT --- Estimate parameters for each value of CELLINE_DRUG_ADDITIVE_ADDITIVECONC (GROUPER, previous line) = combination of the covariates, with REPEAT being the IIV and REP in the RUV
       DROP ; CELLINE_NAME
       DROP ; DRUG_NAME
       DROP ; ADDITIVE_NAME
       DROP ; CELLINE_DRUG_ADDITIVE_NAME

$DATA onco_combine_therapy.csv IGNORE=@

$PRED ;;; 2 CELL LINEs
      	IF(CELLINE .EQ. 1)        TOP_CELLINE = THETA(1)  ;reference --- 1 FIX
      	IF(CELLINE .EQ. 1)     BOTTOM_CELLINE = THETA(2)  ;reference --- 1 FIX
      	IF(CELLINE .EQ. 1)       IC50_CELLINE = THETA(3)  ;reference --- 1 FIX
      	IF(CELLINE .EQ. 1)      GAMMA_CELLINE = THETA(4)  ;reference --- 1 FIX

      	IF(CELLINE .EQ. 2)        TOP_CELLINE = THETA(5)
      	IF(CELLINE .EQ. 2)     BOTTOM_CELLINE = THETA(6)
      	IF(CELLINE .EQ. 2)       IC50_CELLINE = THETA(7)
      	IF(CELLINE .EQ. 2)      GAMMA_CELLINE = THETA(8)

      
      ;;; 2 active drugs
      	IF(DRUG .EQ. 1)    TOP_DRUG = THETA(9) 		;reference --- 1 FIX
      	IF(DRUG .EQ. 1) BOTTOM_DRUG = THETA(10) 	;reference --- 1 FIX
      	IF(DRUG .EQ. 1)   IC50_DRUG = THETA(11) 	;reference --- 1 FIX
      	IF(DRUG .EQ. 1)  GAMMA_DRUG = THETA(12) 	;reference --- 1 FIX

	IF(DRUG .EQ. 2)    TOP_DRUG = THETA(13) 	;reference --- 1 FIX
      	IF(DRUG .EQ. 2) BOTTOM_DRUG = THETA(14) 	;reference --- 1 FIX
      	IF(DRUG .EQ. 2)   IC50_DRUG = THETA(15) 	;reference --- 1 FIX
      	IF(DRUG .EQ. 2)  GAMMA_DRUG = THETA(16) 	;reference --- 1 FIX


      ;;; ADDITIVE ATUX effect
      	IF(ADDITIVECONC .EQ. 0)    TOP_ADDITIVE = THETA(17)	; no ADDITIVE/HAP/SMAP is the same, can be pooled!
      	IF(ADDITIVECONC .EQ. 0) BOTTOM_ADDITIVE = THETA(18)	; no ADDITIVE/HAP/SMAP is the same, can be pooled!
      	IF(ADDITIVECONC .EQ. 0)   IC50_ADDITIVE = THETA(19)	; no ADDITIVE/HAP/SMAP is the same, can be pooled!
      	IF(ADDITIVECONC .EQ. 0)  GAMMA_ADDITIVE = THETA(20)	; no ADDITIVE/HAP/SMAP is the same, can be pooled!

	IF(ADDITIVECONC .EQ. 3)    TOP_ADDITIVE = THETA(21)	
      	IF(ADDITIVECONC .EQ. 3) BOTTOM_ADDITIVE = THETA(22)	
      	IF(ADDITIVECONC .EQ. 3)   IC50_ADDITIVE = THETA(23)	
      	IF(ADDITIVECONC .EQ. 3)  GAMMA_ADDITIVE = THETA(24)

	IF(ADDITIVECONC .EQ. 4.5)    TOP_ADDITIVE = THETA(25)	
      	IF(ADDITIVECONC .EQ. 4.5) BOTTOM_ADDITIVE = THETA(26)	
      	IF(ADDITIVECONC .EQ. 4.5)   IC50_ADDITIVE = THETA(27)	
      	IF(ADDITIVECONC .EQ. 4.5)  GAMMA_ADDITIVE = THETA(28)      

	IF(ADDITIVECONC .EQ. 6)    TOP_ADDITIVE = THETA(29)
      	IF(ADDITIVECONC .EQ. 6)    BOTTOM_ADDITIVE = THETA(30)
      	IF(ADDITIVECONC .EQ. 6)    IC50_ADDITIVE = THETA(31)
      	IF(ADDITIVECONC .EQ. 6)    GAMMA_ADDITIVE = THETA(32)


	;;; inter-day variability
	DAY1=0
	DAY2=0
	DAY3=0
	DAY4=0


	IF(DAY.EQ.1) DAY1=1
	IF(DAY.EQ.2) DAY2=1
	IF(DAY.EQ.3) DAY3=1
	IF(DAY.EQ.4) DAY4=1


	TVTOP    =    TOP_CELLINE *    TOP_DRUG *    TOP_ADDITIVE
 	 TOP    =    TVTOP * EXP(ETA(1)   + DAY1 * ETA(5)  + DAY2 * ETA(6)  + DAY3 * ETA(7) + DAY4 * ETA(8))  

	TVBOTTOM = BOTTOM_CELLINE * BOTTOM_DRUG * BOTTOM_ADDITIVE
  	BOTTOM = TVBOTTOM * EXP(ETA(2)    + DAY1 * ETA(9) + DAY2 * ETA(10) + DAY3 * ETA(11) + DAY4 * ETA(12)) 

	TVIC50   =   IC50_CELLINE *   IC50_DRUG *   IC50_ADDITIVE
  	IC50   =   TVIC50 * EXP(ETA(3)    + DAY1 * ETA(13) + DAY2 * ETA(14) + DAY3 * ETA(15) + DAY4 * ETA(16)) 

	TVGAMMA  =  GAMMA_CELLINE *  GAMMA_DRUG *  GAMMA_ADDITIVE
  	GAMMA  =  TVGAMMA * EXP(ETA(4) + DAY1 * ETA(17) + DAY2 * ETA(18) + DAY3 * ETA(19) + DAY4 * ETA(20)) 

	IMAX = (TOP-BOTTOM)/TOP

	IPRED = BOTTOM + (TOP-BOTTOM) / (1 + (((DRUGCONC) ** GAMMA)/(IC50 ** GAMMA))) + 0.00001

	W = SQRT(IPRED ** 2 * SIGMA(1,1) + SIGMA(2,2))
	Y = IPRED + (IPRED * EPS(1)) + EPS(2)

	IRES = DV - IPRED
	IWRES = IRES / W

$THETA
(1) FIX ;reference cell line WT THETA(1)
(1) FIX ;reference cell line WT
(1) FIX ;reference cell line WT
(1) FIX ;reference cell line WT
(0, 1.02) ;THETA(5)
(0, 1.19) ;
(0, 1.17) ;
(0, 1.03) ;
(1) FIX ;reference Dactolisib THETA(9)
(1) FIX ;reference Dactolisib
(1) FIX ;reference Dactolisib
(1) FIX ;reference Dactolisib
(0, 1.24) ;THETA(13)
(0, 2.34) ;
(0, 206) ;
(0, 1.38) ;
(0, 0.987) ;THETA(17)
(0, 0.0978) ;
(0, 0.0166) ;
(0, 1.07) ;
(0, 0.887) ;THETA(21)
(0, 0.0898) ;
(0, 0.015) ;
(0, 0.96) ;
(0, 0.948) ;THETA(25)
(0, 0.106) ;
(0, 0.016) ;
(0, 1.03) ;
(0, 0.79) ;THETA(29)
(0, 0.11) ;
(0, 0.0191) ;
(0, 1) ;

$OMEGA	;;; inter-experiment variability
 0 FIX 
 4.27E-07
 0.179
 5.15E-07 ; 
$OMEGA BLOCK(1) 0.0388 ; inter-day variability TOP  = intra-experiment variability
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) 0.114 ; inter-day variability BOTTOM
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) 0.0942 ; inter-day variability IC50
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) 5.08E-06 ; inter-day variability GAMMA
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME

$SIGMA ;;; residual variability
 0.0329
 0 FIX 

$COVARIANCE PRINT=E UNCONDITIONAL ;MATRIX=R

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=999999 POSTHOC PRINT=5 NOABORT SADDLE_RESET=1 ;NSIG=3 SIGL=9

$TABLE ID CELLINE ADDITIVE ADDITIVECONC DRUG DAY REP DRUGCONC
DV Y PRED IPRED IRES IWRES CWRES
TOP BOTTOM IC50 GAMMA
TVTOP TVBOTTOM TVIC50 TVGAMMA
IMAX ETA(1) ETA(2) ETA(3) ETA(4)
ONEHEADER NOPRINT FILE = onco_combine_therapy_out_002.csv

