//undefine the grid size parameters if they have already been defined
#ifdef N1 
#undef N1
#endif

#ifdef N2 
#undef N2
#endif

#ifdef N3 
#undef N3
#endif

#define N1 64
#define N2 64
#define N3 1

// atch adjusts
#undef WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS
#undef WENO_USE_PRIM_REDUCTION
#undef LIMIT_FLUXC2A_PRIM_CHANGE
#undef COMPDIM
#undef MAXBND
#undef MCOORD
#undef PRODUCTION
#undef FULLOUTPUT
#undef MAILWHENDONE
#undef EMAILMESSAGE
#undef EMAILADDRESS 
#undef PERFTEST
#undef DOAVG
#undef DOJETDIAG
#undef DOAVG2
#undef DODEBUG
#undef DO_WENO_DEBUG
#undef DOENODEBUG
#undef DODISS
#undef DOLUMVSR
#undef DODISSVSR
#undef DOFIELDLINE
#undef ROEAVERAGEDWAVESPEED
#undef ATHENAROE
#undef STOREWAVESPEEDS
#undef USESTOREDSPEEDSFORFLUX
#undef VCHARTYPE
#undef PRECISEINVERSION
#undef WHICHVEL
#undef WHICHEOM
#undef REMOVERESTMASSFROMUU
#undef RELTYPE
#undef EOMTYPE
#undef UTOPRIMTRYAGAIN
#undef DOENTROPY
#undef WHICHENTROPYEVOLVE
#undef FIXUPAFTERINIT
#undef FIXUPAFTERRESTART
#undef CHECKSOLUTION
#undef GAMMAPERCDIFFMAX
#undef UPERCDIFFMAX
#undef DOEXTRAINTERP
#undef LIMADJUST
#undef HYDROLIMADJUSTONLY
#undef FLUXADJUST
#undef HYDROFLUXADJUSTONLY
#undef STEPOVERNEGU
#undef STEPOVERNEGRHO
#undef STEPOVERNEGRHOU
#undef UTOPRIMADJUST
#undef UTOPRIMFAILRETURNTYPE
#undef COORDSINGFIX
#undef SINGSMALL
#undef VOLUMEDIFF
#undef MINDT
#undef JONCHECKS
#undef JONCHECKS2
#undef FLOORDIAGS
#undef ANALYTICCONNECTION
#undef ANALYTICSOURCE
#undef OUTFLOWAVOIDBC
#undef FLUXDIMENSPLIT
#undef A2CDIMENSPLIT
#undef DODQMEMORY
#undef DOENOFLUXMEMORY
#undef BOUNDARYINTERPADJUST
#undef COMPUTEFRDOT
#undef CALCFARADAYANDCURRENTS
#undef WHICHCURRENTCALC
#undef FARADAYT0
#undef CURRENTST0


#undef EVOLVECHECKS
#undef FIXUPZONES
#undef HLLBOUNDARY
#undef FIXUPFLUX
#undef ZEROOUTFLOWFLUX
#undef ZEROPOLEFLUX
#undef RESCALEINTERP
#undef BDIRCONT
#undef HYPERHLL
#undef HORIZONSUPERFAST
#undef VARTOINTERP
#undef VARTOINTERPFIELD



#undef USEAVGPRIMITIVEFORWENOFLAT
#undef USEPRIMITIVEFROMAVGCONSERVED
#undef CONTACTINDICATOR
#undef COMPUTEDRHODP
#undef SUPERFASTDIVREDUCE
#undef MINPREFORDER
#undef SHOCKINDICATOR
#undef WHICHPARA
#undef HOT2COLD



// jon adjusts
#define WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS 1
//#define WENO_USE_LIM_PRIM_CORRECTION_FOR_FLUX_LIMITING 1

#define WENO_USE_PRIM_REDUCTION 1

#define LIMIT_FLUXC2A_PRIM_CHANGE 0



#define MCOORD HTMETRICACCURATE
//#define MCOORD SPCMINKMETRIC
#define COMPDIM 3
#define MAXBND 3
#define PRODUCTION 0
#define FULLOUTPUT MAXBND

#define MAILWHENDONE 1
#define MAILFROMREMOTE 0
#define REMOTEHOST "bh.astro.uiuc.edu"
#define EMAILADDRESS "jmckinney@cfa.harvard.edu"
#define EMAILMESSAGE "Done with GRMHD run DEFAULT"
#define PERFTEST 0
#define DOAVG 0
#define DOJETDIAG 1
#define DOAVG2 0
#define DODEBUG 1
#define DO_WENO_DEBUG 0
#define DOENODEBUG 0
#define DODISS 0
#define DOLUMVSR 1
#define DODISSVSR 0
#define DOFIELDLINE 1
#define ROEAVERAGEDWAVESPEED 0
#define ATHENAROE 0
#define STOREWAVESPEEDS 0  //set this and the following one to unity to use the DONOR interpolated states for computing wavespeeds
#define USESTOREDSPEEDSFORFLUX 0
#define VCHARTYPE VERYLOCALVCHAR
#define PRECISEINVERSION 1
#define WHICHVEL VELREL4
#define WHICHEOM WITHGDET
#define REMOVERESTMASSFROMUU 2
#define RELTYPE RELEOM
//#define EOMTYPE EOMGRMHD
//#define EOMTYPE EOMFFDE
#define EOMTYPE EOMCOLDGRMHD
#define UTOPRIMTRYAGAIN 0
#define DOENTROPY DONOENTROPY
#define WHICHENTROPYEVOLVE EVOLVESIMPLEENTROPY
#define FIXUPAFTERINIT 1
#define FIXUPAFTERRESTART 1
#define CHECKSOLUTION 0
#define GAMMAPERCDIFFMAX 2.0
#define UPERCDIFFMAX 10.0
#define LIMADJUST LIMITERFIXED
#define HYDROLIMADJUSTONLY 0
#define FLUXADJUST FLUXFIXED
#define HYDROFLUXADJUSTONLY 0
#define STEPOVERNEGU -1
#define STEPOVERNEGRHO -1
#define STEPOVERNEGRHOU -1
#define UTOPRIMADJUST UTOPRIMAVG
#define UTOPRIMFAILRETURNTYPE UTOPRIMRETURNADJUSTED
#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to
//#define SINGSMALL (1E-3)
#define SINGSMALL (1E-14) // must be larger than machine precision to work for outer M_PI boundary!
// Hawley uses 0.06283 (0.02Pi)

#define VOLUMEDIFF 0
#define MINDT 1.e-20 
#define JONCHECKS 1    //SASMARK - do I need this?
#define JONCHECKS2 1   //SASMARK - do I need this?
#define FLOORDIAGS 1
#define ANALYTICCONNECTION 0  //SASMARK - Don't I need this?
#define ANALYTICSOURCE 0
#define OUTFLOWAVOIDBC 0
#define FLUXDIMENSPLIT QUASISTRANG
#define A2CDIMENSPLIT QUASISTRANG
#define DODQMEMORY 1
#define DOENOFLUXMEMORY 1
#define BOUNDARYINTERPADJUST 0  //should be set to zero always
#define COMPUTEFRDOT 0
#define CALCFARADAYANDCURRENTS 1
#define WHICHCURRENTCALC CURRENTCALC1
#define FARADAYT0 1
#define CURRENTST0 1


#define EVOLVECHECKS 0
#define FIXUPZONES FIXUP1ZONE
#define HLLBOUNDARY 0
#define FIXUPFLUX 0
#define ZEROOUTFLOWFLUX 0
#define ZEROPOLEFLUX 0
#define BDIRCONT 1
#define HYPERHLL 0
#define HORIZONSUPERFAST 0

#define VARTOINTERP PRIMTOINTERP
//#define VARTOINTERP PRIMTOINTERP_RHOU
//#define VARTOINTERP PRIMTOINTERP_VSQ
#define RESCALEINTERP 1
#define DOEXTRAINTERP 0

#define NOSPECIALFIELD 0
#define PULSARFIELD 1

#define VARTOINTERPFIELD PULSARFIELD

#define USEAVGPRIMITIVEFORWENOFLAT 1
#define USEPRIMITIVEFROMAVGCONSERVED 0
#define CONTACTINDICATOR 0
#define COMPUTEDRHODP 1
#define SUPERFASTDIVREDUCE 0
#define MINPREFORDER 3
#define SHOCKINDICATOR 1
#define WHICHPARA PARA4

#undef DO_VORTICITY_IMAGE
#define DO_VORTICITY_IMAGE 0

#define HOT2COLD 0


struct Ccoordparams {
  double timescalefactor;
}  coordparams;
