;; x1. Author: apple

;; 3. Label:

$SIZES LTH=9999 LVR=9999
       DIMCNS=9999 DIMTMP=9999 DIMNEW=9999
;        DIMVRB=9999 LIM1=9999 MAXFCN=9999 NO=9999 MAXPTHETA=9999
$PROBLEM Quantitative drug sensitivity scoring of 1 coumpound in 5 strains

$INPUT CELLINE         ; 5 cell lines: EGFR, WT, R183W, P179R, S256F
       DRUG            ; Clofarabine
       ADDITIVE        ; OA
       ADDITIVECONC    ; 0,10 nM
       DAY             ; inter-day repeats
       REP             ; within-day repeats
       DRUGCONC        ; 0, 0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 
       DV           
       DROP ; ID option 1 --- CELLINE_DRUG_ADDITIVE
       ID   ; ID option 2 --- CELLINE_DRUG_ADDITIVE_ADDITIVECONC ; --- Estimate parameters for each value of CELLINE_DRUG_ADDITIVE_ADDITIVECONC = combination of the covariates, with REPEAT in the IOV/RUV and REP in the RUV, or ... define CELLINE_DRUG_ADDITIVE_ADDITIVECONC_DAY as ID ...
       DROP ; ID option 3 --- CELLINE_DRUG_ADDITIVE_ADDITIVECONC_DAY ; CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT --- Estimate parameters for each value of CELLINE_DRUG_ADDITIVE_ADDITIVECONC (GROUPER, previous line) = combination of the covariates, with REPEAT being the IIV and REP in the RUV
       DROP ; CELLINE_NAME
       DROP ; DRUG_NAME
       DROP ; ADDITIVE_NAME
       DROP ; CELLINE_DRUG_ADDITIVE_NAME

$DATA DS_2023_final.csv IGNORE=@

$PRED ;;; CELL LINE
      IF(CELLINE .EQ. 1)        TOP_CELLINE = THETA(1)  ;reference --- 1 FIX
      IF(CELLINE .EQ. 1)     BOTTOM_CELLINE = THETA(2)  ;reference --- 1 FIX
      IF(CELLINE .EQ. 1)       IC50_CELLINE = THETA(3)  ;reference --- 1 FIX
      IF(CELLINE .EQ. 1)      GAMMA_CELLINE = THETA(4)  ;reference --- 1 FIX

      IF(CELLINE .EQ. 2)        TOP_CELLINE = THETA(5)
      IF(CELLINE .EQ. 2)     BOTTOM_CELLINE = THETA(6)
      IF(CELLINE .EQ. 2)       IC50_CELLINE = THETA(7)
      IF(CELLINE .EQ. 2)      GAMMA_CELLINE = THETA(8)

      IF(CELLINE .EQ. 3)        TOP_CELLINE = THETA(9)
      IF(CELLINE .EQ. 3)     BOTTOM_CELLINE = THETA(10)
      IF(CELLINE .EQ. 3)       IC50_CELLINE = THETA(11)
      IF(CELLINE .EQ. 3)      GAMMA_CELLINE = THETA(12)
	  
	  IF(CELLINE .EQ. 4)        TOP_CELLINE = THETA(13)
      IF(CELLINE .EQ. 4)     BOTTOM_CELLINE = THETA(14)
      IF(CELLINE .EQ. 4)       IC50_CELLINE = THETA(15)
      IF(CELLINE .EQ. 4)      GAMMA_CELLINE = THETA(16)

      IF(CELLINE .EQ. 5)        TOP_CELLINE = THETA(17)
      IF(CELLINE .EQ. 5)     BOTTOM_CELLINE = THETA(18)
      IF(CELLINE .EQ. 5)       IC50_CELLINE = THETA(19)
      IF(CELLINE .EQ. 5)      GAMMA_CELLINE = THETA(20)

      ;;; 1 base drug
      IF(DRUG .EQ. 1)    TOP_DRUG = THETA(21) ;reference --- 1 FIX
      IF(DRUG .EQ. 1) BOTTOM_DRUG = THETA(22) ;reference --- 1 FIX
      IF(DRUG .EQ. 1)   IC50_DRUG = THETA(23) ;reference --- 1 FIX
      IF(DRUG .EQ. 1)  GAMMA_DRUG = THETA(24) ;reference --- 1 FIX

      ;;; ADDITIVE OA effect
      IF(ADDITIVECONC .EQ. 0)    TOP_ADDITIVE = THETA(25) ; no ADDITIVE/HAP/SMAP is the same, can be pooled!
      IF(ADDITIVECONC .EQ. 0) BOTTOM_ADDITIVE = THETA(26) ; no ADDITIVE/HAP/SMAP is the same, can be pooled!
      IF(ADDITIVECONC .EQ. 0)   IC50_ADDITIVE = THETA(27) ; no ADDITIVE/HAP/SMAP is the same, can be pooled!
      IF(ADDITIVECONC .EQ. 0)  GAMMA_ADDITIVE = THETA(28) ; no ADDITIVE/HAP/SMAP is the same, can be pooled!

      IF(ADDITIVECONC .EQ. 10)    TOP_ADDITIVE = THETA(29)
      IF(ADDITIVECONC .EQ. 10)    BOTTOM_ADDITIVE = THETA(30)
      IF(ADDITIVECONC .EQ. 10)    IC50_ADDITIVE = THETA(31)
      IF(ADDITIVECONC .EQ. 10)    GAMMA_ADDITIVE = THETA(32)

  

; inter-day variability
DAY1=0
DAY2=0
DAY3=0
DAY4=0
DAY5=0
DAY6=0
DAY7=0

IF(DAY.EQ.1) DAY1=1
IF(DAY.EQ.2) DAY2=1
IF(DAY.EQ.3) DAY3=1
IF(DAY.EQ.4) DAY4=1
IF(DAY.EQ.5) DAY5=1
IF(DAY.EQ.6) DAY6=1
IF(DAY.EQ.7) DAY7=1



TVTOP    =    TOP_CELLINE *    TOP_DRUG *    TOP_ADDITIVE
  TOP    =    TVTOP *EXP(ETA(1)    + DAY1*ETA(5)  + DAY2*ETA(6)  + DAY3*ETA(7))  

TVBOTTOM = BOTTOM_CELLINE * BOTTOM_DRUG * BOTTOM_ADDITIVE
  BOTTOM = TVBOTTOM *EXP(ETA(2)    + DAY1*ETA(8) + DAY2*ETA(9) + DAY3*ETA(10)) 

TVIC50   =   IC50_CELLINE *   IC50_DRUG *   IC50_ADDITIVE
  IC50   =   TVIC50 *EXP(ETA(3)    + DAY1*ETA(11) + DAY2*ETA(12) + DAY3*ETA(13)) 

TVGAMMA  =  GAMMA_CELLINE *  GAMMA_DRUG *  GAMMA_ADDITIVE
  GAMMA  =  TVGAMMA *EXP(ETA(4) + DAY1*ETA(14) + DAY2*ETA(15) + DAY3*ETA(16)) 

IMAX = (TOP-BOTTOM)/TOP

IPRED = BOTTOM + (TOP-BOTTOM)/(1+(((DRUGCONC)**GAMMA)/(IC50**GAMMA)))+0.00001

W = SQRT(IPRED**2*SIGMA(1,1) + SIGMA(2,2))
Y = IPRED + (IPRED*EPS(1)) + EPS(2)

IRES = DV-IPRED
IWRES = IRES/W

$THETA
(1) FIX ; reference cell line EFGR (1)
(1) FIX ; reference cell line EFGR
(1) FIX ; reference cell line EFGR
(1) FIX ; reference cell line EFGR
(0, 1.01) ;5
(0, 1.13) ;
(0, 1.11) ;
(0, 1.01) ;
(0, 0.969) ;9
(0, 0.995) ;
(0, 1.11) ;
(0, 1.09) ;
(0, 1.01) ;13
(0, 1.13) ;
(0, 1.11) ;
(0, 1.01) ;
(0, 0.969) ;17
(0, 0.995) ;
(0, 1.11) ;
(0, 1.09) ;
(1) FIX ; reference Clofarabine (21)
(1) FIX ; reference Clofarabine
(1) FIX ; reference Clofarabine
(1) FIX ; reference Clofarabine
(0, 1.04) ;25
(0, 0.409) ;
(0, 2.26) ;
(0, 0.647) ;
(0, 1.03) ;29
(0, 0.449) ;
(0, 0.14) ;
(0, 0.969) ;

$OMEGA       ; inter-experiment variability
 0 FIX 
 0.00427
 0.0355
 0.00515 ; 
$OMEGA BLOCK(1) 0.00912 ; inter-day variability TOP  = intra-experiment variability
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) 0.0542 ; inter-day variability BOTTOM
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) 0.0524 ; inter-day variability IC50
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) 0.0508 ; inter-day variability GAMMA
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME
$OMEGA BLOCK(1) SAME

$SIGMA 0.00986 ; residual variability
 0 FIX 

$COVARIANCE PRINT=E UNCONDITIONAL ;MATRIX=R

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=999999 POSTHOC PRINT=5 NOABORT SADDLE_RESET=1 ;NSIG=3 SIGL=9

$TABLE ID CELLINE ADDITIVE ADDITIVECONC DRUG DAY REP DRUGCONC
DV Y PRED IPRED IRES IWRES CWRES
TOP BOTTOM IC50 GAMMA
TVTOP TVBOTTOM TVIC50 TVGAMMA
IMAX ETA(1) ETA(2) ETA(3) ETA(4)
ONEHEADER NOPRINT FILE=runfinal_2023_out.csv
