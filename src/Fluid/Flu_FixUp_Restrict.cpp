#include "GAMER.h"




#if ( MODEL == ELBDM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Restrict4_Weights
// Description :  1D weights W[i][s] for the fourth-order restriction of PS2 fine cells onto PS1 coarse cells (point values at the
//                coarse cell centres, i.e. in between the fine cells 2i and 2i+1)
//                --> interior: cubic Lagrange interpolation (-1, 9, 9, -1)/16 on the fine cells 2i-1 ... 2i+2
//                --> first and last coarse cell: one-sided cubic stencils (5, 15, -5, 1)/16 and (1, -5, 15, 5)/16
//-------------------------------------------------------------------------------------------------------
static void Restrict4_Weights( real W[PS1][PS2] )
{
   for (int i=0; i<PS1; i++)
   for (int s=0; s<PS2; s++)   W[i][s] = (real)0.0;

   for (int i=1; i<PS1-1; i++)
   {
      W[i][2*i-1] = (real)(-1.0/16.0);
      W[i][2*i  ] = (real)( 9.0/16.0);
      W[i][2*i+1] = (real)( 9.0/16.0);
      W[i][2*i+2] = (real)(-1.0/16.0);
   }

   W[0][0] = (real)( 5.0/16.0);   W[0][1] = (real)(15.0/16.0);   W[0][2] = (real)(-5.0/16.0);   W[0][3] = (real)( 1.0/16.0);

   W[PS1-1][PS2-4] = (real)( 1.0/16.0);   W[PS1-1][PS2-3] = (real)(-5.0/16.0);
   W[PS1-1][PS2-2] = (real)(15.0/16.0);   W[PS1-1][PS2-1] = (real)( 5.0/16.0);

} // FUNCTION : Restrict4_Weights



//-------------------------------------------------------------------------------------------------------
// Function    :  Restrict4_Block
// Description :  Separable fourth-order restriction of a PS2^3 block of fine cells onto the PS1^3 cells of the father patch
//-------------------------------------------------------------------------------------------------------
static void Restrict4_Block( real Dst[PS1][PS1][PS1], const real Src[PS2][PS2][PS2], const real W[PS1][PS2] )
{
   real T1[PS2][PS2][PS1], T2[PS2][PS1][PS1];

// x direction
   for (int k=0; k<PS2; k++)
   for (int j=0; j<PS2; j++)
   for (int i=0; i<PS1; i++)
   {
      real Sum = (real)0.0;
      for (int s=0; s<PS2; s++)   if ( W[i][s] != (real)0.0 )   Sum += W[i][s]*Src[k][j][s];
      T1[k][j][i] = Sum;
   }

// y direction
   for (int k=0; k<PS2; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
      real Sum = (real)0.0;
      for (int s=0; s<PS2; s++)   if ( W[j][s] != (real)0.0 )   Sum += W[j][s]*T1[k][s][i];
      T2[k][j][i] = Sum;
   }

// z direction
   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
      real Sum = (real)0.0;
      for (int s=0; s<PS2; s++)   if ( W[k][s] != (real)0.0 )   Sum += W[k][s]*T2[s][j][i];
      Dst[k][j][i] = Sum;
   }

} // FUNCTION : Restrict4_Block



//-------------------------------------------------------------------------------------------------------
// Function    :  Restrict4_Gather
// Description :  Copy variable "v" of the eight son patches SonPID0 ... SonPID0+7 into one PS2^3 block
//-------------------------------------------------------------------------------------------------------
static void Restrict4_Gather( real Blk[PS2][PS2][PS2], const int SonLv, const int SonPID0, const int SonFluSg, const int v )
{
   for (int LocalID=0; LocalID<8; LocalID++)
   {
      const int Ox = TABLE_02( LocalID, 'x', 0, PS1 );
      const int Oy = TABLE_02( LocalID, 'y', 0, PS1 );
      const int Oz = TABLE_02( LocalID, 'z', 0, PS1 );

      const real (*Son)[PS1][PS1] = amr->patch[SonFluSg][SonLv][SonPID0+LocalID]->fluid[v];

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)   Blk[Oz+k][Oy+j][Ox+i] = Son[k][j][i];
   }

} // FUNCTION : Restrict4_Gather



//-------------------------------------------------------------------------------------------------------
// Function    :  Restrict4_BoxAverage
// Description :  Box average (second-order restriction) of cell (i,j,k) of the father patch, used as positivity fallback
//-------------------------------------------------------------------------------------------------------
static real Restrict4_BoxAverage( const real Src[PS2][PS2][PS2], const int i, const int j, const int k )
{
   real Sum = (real)0.0;

   for (int dk=0; dk<2; dk++)
   for (int dj=0; dj<2; dj++)
   for (int di=0; di<2; di++)   Sum += Src[2*k+dk][2*j+dj][2*i+di];

   return (real)0.125*Sum;

} // FUNCTION : Restrict4_BoxAverage
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_FixUp_Restrict
// Description :  Replace the data at level "FaLv" by the average data at level "FaLv+1"
//
// Note        :  1. Use the input parameters "TVarCC" and "TVarFC" to determine the target cell- and
//                   face-centered variables on level "FaLv", respectively
//                2. For MHD, this routine currently always restrict all three B field components
//                   --> Do not distinguish _MAGX, _MAGY, _MAGZ, and _MAG in TVarFC
//                3. Invoked by EvolveLevel()
//                4. ELBDM_HYBRID + LOAD_BALANCE: Backward matching of phase field for ELBDM_MATCH_PHASE requires OPT__LB_EXCHANGE_FATHER
//
// Parameter   :  FaLv     : Target refinement level at which the data are going to be replaced
//                SonFluSg : Fluid sandglass at level "FaLv+1"
//                FaFluSg  : Fluid sandglass at level "FaLv"
//                SonMagSg : B field sandglass at level "FaLv+1"
//                FaMagSg  : B field sandglass at level "FaLv"
//                SonPotSg : Potential sandglass at level "FaLv+1"
//                FaPotSg  : Potential sandglass at level "FaLv"
//                TVarCC   : Target cell-centered variables on level "FaLv"
//                           --> Supported variables in different models:
//                               HYDRO        : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY [, _POTE] [, BIDX(field_index)]
//                               ELBDM_WAVE   : _DENS, _REAL, _IMAG, [, _POTE]
//                               ELBDM_HYBRID : _DENS, _PHAS [, _POTE]
//                           --> _FLUID, _PASSIVE, and _TOTAL apply to all models
//                TVarFC   : Target face-centered variables on level "FaLv"
//                            --> Supported variables in different models:
//                                HYDRO with MHD : _MAGX, _MAGY, _MAGZ, _MAG
//                                ELBDM          : none
//                            --> But it currently does not distinguish _MAGX, _MAGY, _MAGZ, and _MAG
//-------------------------------------------------------------------------------------------------------
void Flu_FixUp_Restrict( const int FaLv, const int SonFluSg, const int FaFluSg, const int SonMagSg, const int FaMagSg,
                         const int SonPotSg, const int FaPotSg, const long TVarCC, const long TVarFC )
{

   const int SonLv = FaLv + 1;


// check
   if ( FaLv < 0  ||  FaLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );

   if (  ( TVarCC & _TOTAL )  &&  ( SonFluSg != 0 && SonFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonFluSg", SonFluSg );

   if (  ( TVarCC & _TOTAL )  &&  ( FaFluSg != 0 && FaFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaFluSg", FaFluSg );

#  ifdef GRAVITY
   if (  ( TVarCC & _POTE )  &&  ( SonPotSg != 0 && SonPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonPotSg", SonPotSg );

   if (  ( TVarCC & _POTE )  &&  ( FaPotSg != 0 && FaPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaPotSg", FaPotSg );
#  endif

#  ifdef MHD
   if (  ( TVarFC & _MAG )  &&  TVarFC != _MAG  )
      Aux_Error( ERROR_INFO, "must work on all three magnetic components at once !!\n" );

   if (  ( TVarFC & _MAG )  &&  ( SonMagSg != 0 && SonMagSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonMagSg", SonMagSg );

   if (  ( TVarFC & _MAG )  &&  ( FaMagSg != 0 && FaMagSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaMagSg", FaMagSg );
#  endif


// nothing to do on the top level
   if ( FaLv == TOP_LEVEL )
   {
      Aux_Message( stderr, "WARNING : applying \"%s\" to the top level is meaningless !!\n", __FUNCTION__ );
      return;
   }

// nothing to do if there is no real patch at lv+1
   if ( amr->NPatchComma[SonLv][1] == 0 )    return;

// check the synchronization
   Mis_CompareRealValue( Time[FaLv], Time[SonLv], __FUNCTION__, true );


   const bool ResFlu  = TVarCC & _TOTAL;
#  ifdef GRAVITY
   const bool ResPot  = TVarCC & _POTE;
#  else
   const bool ResPot  = false;
#  endif
#  ifdef MHD
   const bool ResMag  = TVarFC & _MAG;
#  else
   const bool ResMag  = false;
#  endif
   const int PS1_half = PS1 / 2;


// determine the indices of fluid components to be restricted: TFluVarIdx = [0 ... NCOMP_TOTAL-1]
   int NFluVar=0, TFluVarIdxList[NCOMP_TOTAL];

#  if ( MODEL == ELBDM )
// set correct restriction options for wave-wave, fluid-wave and fluid-fluid level phase restriction
// wave-wave:   phase restriction on if OPT__RES_PHASE is set
// fluid-wave:  always use phase restriction
// fluid-fluid: treat phase restriction the same as normal restriction
// --> ELBDM_MATCH_PHASE is only applied to fluid-wave
   const bool ResWave  = ( TVarCC & _REAL ) || ( TVarCC & _IMAG );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
// restrict phase during wave-wave-level restriction if OPT__RES_PHAS is enabled
   const bool ResWWPha = OPT__RES_PHASE && ResWave &&  amr->use_wave_flag[FaLv] && amr->use_wave_flag[SonLv];
// always restrict phase during fluid-wave-level restriction
   const bool ResWFPha =                   ResWave && !amr->use_wave_flag[FaLv] && amr->use_wave_flag[SonLv];
   const bool ResPha   = ResWWPha || ResWFPha;
#  else
// restrict phase if OPT__RES_PHAS is enabled
   const bool ResPha   = OPT__RES_PHASE && ResWave;
#  endif

// update the components to be restricted depending on whether phase restriction is enabled
// --> when applying phase restriction, TFluVarIdxList[] only needs to record components other than density and wave function
   for (int v=0; v<NCOMP_TOTAL; v++) {
      if ( ResPha ) {
         if ( TVarCC & (1L<<v)  &&  v != DENS  &&  v != REAL  &&  v != IMAG )    TFluVarIdxList[ NFluVar ++ ] = v;
      } else {
         if ( TVarCC & (1L<<v)                                              )    TFluVarIdxList[ NFluVar ++ ] = v;
      }
   }

#  else // #if ( MODEL == ELBDM )

   const bool ResPha = false;

   for (int v=0; v<NCOMP_TOTAL; v++)
      if ( TVarCC & (1L<<v) )    TFluVarIdxList[ NFluVar ++ ] = v;
#  endif // #if ( MODEL == ELBDM ) ... else ...


// return immediately if no target variable is found
   if ( NFluVar == 0  &&  !ResPot  &&  !ResMag  &&  !ResPha )
   {
      Aux_Message( stderr, "WARNING : no target variable is found in %s !!\n", __FUNCTION__ );
      return;
   }


// fourth-order restriction (ELBDM only)
#  if ( MODEL == ELBDM )
   const bool Res4     = ( OPT__RESTRICT_ORDER >= 3 );                                  // 3: 4th order for phase/wave function only, 4: all, 5: density only
   const bool Res4Dens = ( OPT__RESTRICT_ORDER == 4  ||  OPT__RESTRICT_ORDER == 5 );
   const bool Res4Pha  = ( OPT__RESTRICT_ORDER == 3  ||  OPT__RESTRICT_ORDER == 4 );
   real Res4W[PS1][PS2];
   if ( Res4 )   Restrict4_Weights( Res4W );
#  else
   const bool Res4     = false;
   const bool Res4Dens = false;
   const bool Res4Pha  = false;
#  endif


// restrict
#  pragma omp parallel for schedule( runtime )
   for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)
   {
      const int FaPID = amr->patch[0][SonLv][SonPID0]->father;

//    check
#     ifdef GAMER_DEBUG
      if ( FaPID < 0 )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d has no father patch (FaPID = %d) !!\n",
                    SonLv, SonPID0, FaPID );

      if ( ResFlu  &&  amr->patch[FaFluSg][FaLv][FaPID]->fluid == NULL )
         Aux_Error( ERROR_INFO, "FaFluSg %d, FaLv %d, FaPID %d has no fluid array allocated !!\n",
                    FaFluSg, FaLv, FaPID );

#     ifdef GRAVITY
      if ( ResPot  &&  amr->patch[FaPotSg][FaLv][FaPID]->pot == NULL )
         Aux_Error( ERROR_INFO, "FaPotSg %d, FaLv %d, FaPID %d has no potential array allocated !!\n",
                    FaPotSg, FaLv, FaPID );
#     endif

#     ifdef MHD
      if ( ResMag  &&  amr->patch[FaMagSg][FaLv][FaPID]->magnetic == NULL )
         Aux_Error( ERROR_INFO, "FaMagSg %d, FaLv %d, FaPID %d has no B field array allocated !!\n",
                    FaMagSg, FaLv, FaPID );
#     endif
#     endif // #ifdef GAMER_DEBUG


//    fourth-order restriction of the fluid data (ELBDM only, OPT__RESTRICT_ORDER = 4)
#     if ( MODEL == ELBDM )
      if ( Res4  &&  ResFlu )
      {
         real (*BlkA)[PS2][PS2] = new real [PS2][PS2][PS2];   // heap: avoid large stack arrays in OpenMP threads
         real (*BlkB)[PS2][PS2] = new real [PS2][PS2][PS2];
         real (*BlkC)[PS2][PS2] = new real [PS2][PS2][PS2];
         real (*DstD)[PS1][PS1] = new real [PS1][PS1][PS1];
         real (*DstP)[PS1][PS1] = new real [PS1][PS1][PS1];

         if ( ResPha )
         {
//          density
            Restrict4_Gather( BlkA, SonLv, SonPID0, SonFluSg, DENS );
            Restrict4_Block ( DstD, BlkA, Res4W );

            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            {
               const real Box = Restrict4_BoxAverage( BlkA, i, j, k );
               if ( !Res4Dens  ||  DstD[k][j][i] <= (real)0.0  ||  DstD[k][j][i] < (real)0.5*Box )   DstD[k][j][i] = Box;
            }

//          phase: unwrap the child phases continuously over the whole block before restricting
            Restrict4_Gather( BlkB, SonLv, SonPID0, SonFluSg, REAL );
            Restrict4_Gather( BlkC, SonLv, SonPID0, SonFluSg, IMAG );

            for (int k=0; k<PS2; k++)
            for (int j=0; j<PS2; j++)
            for (int i=0; i<PS2; i++)
            {
               const real Raw = SATAN2( BlkC[k][j][i], BlkB[k][j][i] );

               if      ( i > 0 )   BlkA[k][j][i] = ELBDM_UnwrapPhase( BlkA[k][j][i-1], Raw );
               else if ( j > 0 )   BlkA[k][j][i] = ELBDM_UnwrapPhase( BlkA[k][j-1][0], Raw );
               else if ( k > 0 )   BlkA[k][j][i] = ELBDM_UnwrapPhase( BlkA[k-1][0][0], Raw );
               else                BlkA[k][j][i] = Raw;
            }

            Restrict4_Block( DstP, BlkA, Res4W );

            if ( !Res4Pha )
            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)   DstP[k][j][i] = Restrict4_BoxAverage( BlkA, i, j, k );

            real (*DFaPtr)[PS1][PS1] = amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS];
            real (*RFaPtr)[PS1][PS1] = amr->patch[FaFluSg][FaLv][FaPID]->fluid[REAL];
            real (*IFaPtr)[PS1][PS1] = amr->patch[FaFluSg][FaLv][FaPID]->fluid[IMAG];
#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            real (*PFaPtr)   [PS1][PS1] = amr->patch[  FaFluSg][FaLv][FaPID]->fluid[PHAS];
            real (*OldPFaPtr)[PS1][PS1] = amr->patch[1-FaFluSg][FaLv][FaPID]->fluid[PHAS];
            if ( amr->FluSgTime[FaLv][1-FaFluSg] < 0.0 )   OldPFaPtr = PFaPtr;
#           endif

            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            {
               const real avgdens = DstD[k][j][i];
               real       avgphase = DstP[k][j][i];

#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               if ( !amr->use_wave_flag[FaLv] && ELBDM_MATCH_PHASE )
                  avgphase = ELBDM_UnwrapPhase( OldPFaPtr[k][j][i], avgphase );

               if ( !amr->use_wave_flag[FaLv] ) {
                  if ( TVarCC & _DENS )   DFaPtr[k][j][i] = avgdens;
                  if ( TVarCC & _PHAS )   PFaPtr[k][j][i] = avgphase;
               } else {
#              endif
                  if ( TVarCC & _DENS )   DFaPtr[k][j][i] = avgdens;
                  if ( TVarCC & _REAL )   RFaPtr[k][j][i] = SQRT(avgdens) * COS(avgphase);
                  if ( TVarCC & _IMAG )   IFaPtr[k][j][i] = SQRT(avgdens) * SIN(avgphase);
#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               } // if ( !amr->use_wave_flag[FaLv] ) ... else ...
#              endif
            }
         } // if ( ResPha )

//       remaining variables
         if ( ResFlu )
         for (int v=0; v<NFluVar; v++)
         {
            const int TFluVarIdx = TFluVarIdxList[v];
            real (*FaPtr)[PS1][PS1] = amr->patch[FaFluSg][FaLv][FaPID]->fluid[TFluVarIdx];

            Restrict4_Gather( BlkA, SonLv, SonPID0, SonFluSg, TFluVarIdx );
            Restrict4_Block ( DstD, BlkA, Res4W );

            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            {
               real Val = DstD[k][j][i];

               const bool IsDens = ( TFluVarIdx == DENS );
               const bool UseBox = IsDens ? !Res4Dens : !Res4Pha;
               const real Box    = Restrict4_BoxAverage( BlkA, i, j, k );

               if ( UseBox )   Val = Box;
               else if ( IsDens  &&  ( Val <= (real)0.0  ||  Val < (real)0.5*Box ) )   Val = Box;

               FaPtr[k][j][i] = Val;
            }
         }

         delete [] BlkA;   delete [] BlkB;   delete [] BlkC;   delete [] DstD;   delete [] DstP;
      } // if ( Res4 )
#     endif // #if ( MODEL == ELBDM )


//    loop over eight sons
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int SonPID = SonPID0 + LocalID;
         const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1_half );
         const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1_half );
         const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1_half );

//       check
#        ifdef GAMER_DEBUG
         if ( ResFlu  &&  amr->patch[SonFluSg][SonLv][SonPID]->fluid == NULL )
            Aux_Error( ERROR_INFO, "SonFluSg %d, SonLv %d, SonPID %d has no fluid array allocated !!\n",
                       SonFluSg, SonLv, SonPID );

#        ifdef GRAVITY
         if ( ResPot  &&  amr->patch[SonPotSg][SonLv][SonPID]->pot == NULL )
            Aux_Error( ERROR_INFO, "SonPotSg %d, SonLv %d, SonPID %d has no potential array allocated !!\n",
                       SonPotSg, SonLv, SonPID );
#        endif

#        ifdef MHD
         if ( ResMag  &&  amr->patch[SonMagSg][SonLv][SonPID]->magnetic == NULL )
            Aux_Error( ERROR_INFO, "SonMagSg %d, SonLv %d, SonPID %d has no B field array allocated !!\n",
                       SonMagSg, SonLv, SonPID );
#        endif
#        endif // #ifdef GAMER_DEBUG

#        if ( MODEL == ELBDM )
//       average phase instead of real and imaginary part if ResPha is set
         if ( ResPha  &&  !Res4 ) {

//          D = DENS, R = REAL, I = IMAG, P = PHAS, S = STUB
            const real (*DSonPtr)  [PS1][PS1] = amr->patch[SonFluSg][SonLv][SonPID]->fluid[DENS];
            const real (*RSonPtr)  [PS1][PS1] = amr->patch[SonFluSg][SonLv][SonPID]->fluid[REAL];
            const real (*ISonPtr)  [PS1][PS1] = amr->patch[SonFluSg][SonLv][SonPID]->fluid[IMAG];

                  real (*DFaPtr)   [PS1][PS1] = amr->patch[ FaFluSg][ FaLv][ FaPID]->fluid[DENS];
                  real (*RFaPtr)   [PS1][PS1] = amr->patch[ FaFluSg][ FaLv][ FaPID]->fluid[REAL];
                  real (*IFaPtr)   [PS1][PS1] = amr->patch[ FaFluSg][ FaLv][ FaPID]->fluid[IMAG];

#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
                  real (*PFaPtr)   [PS1][PS1] = amr->patch[  FaFluSg][ FaLv][ FaPID]->fluid[PHAS];
                  real (*OldPFaPtr)[PS1][PS1] = amr->patch[1-FaFluSg][ FaLv][ FaPID]->fluid[PHAS];
//                handle that we do not have data of previous time step during initialisation corresponding to a negative time
                  if ( amr->FluSgTime[FaLv][1-FaFluSg ] < 0.0 ) {
                     OldPFaPtr = PFaPtr;
                  }
#           endif

            int ii, jj, kk, I, J, K, Ip, Jp, Kp;
            real refphase, avgphase, avgdens;

            for (int k=0; k<PS1_half; k++)  {  K = k*2;  Kp = K+1;  kk = k + Disp_k;
            for (int j=0; j<PS1_half; j++)  {  J = j*2;  Jp = J+1;  jj = j + Disp_j;
            for (int i=0; i<PS1_half; i++)  {  I = i*2;  Ip = I+1;  ii = i + Disp_i;

//             take care to match the child phases before averaging
               refphase    = SATAN2  ( ISonPtr[K ][J ][I ], RSonPtr[K ][J ][I ]  );
               avgphase    = 0.125 * (                   refphase                                                    +
                                       ELBDM_UnwrapPhase(refphase, SATAN2(ISonPtr[K ][J ][Ip], RSonPtr[K ][J ][Ip])) +
                                       ELBDM_UnwrapPhase(refphase, SATAN2(ISonPtr[K ][Jp][I ], RSonPtr[K ][Jp][I ])) +
                                       ELBDM_UnwrapPhase(refphase, SATAN2(ISonPtr[Kp][J ][I ], RSonPtr[Kp][J ][I ])) +
                                       ELBDM_UnwrapPhase(refphase, SATAN2(ISonPtr[K ][Jp][Ip], RSonPtr[K ][Jp][Ip])) +
                                       ELBDM_UnwrapPhase(refphase, SATAN2(ISonPtr[Kp][Jp][I ], RSonPtr[Kp][Jp][I ])) +
                                       ELBDM_UnwrapPhase(refphase, SATAN2(ISonPtr[Kp][J ][Ip], RSonPtr[Kp][J ][Ip])) +
                                       ELBDM_UnwrapPhase(refphase, SATAN2(ISonPtr[Kp][Jp][Ip], RSonPtr[Kp][Jp][Ip])) );
               avgdens     = 0.125 * ( DSonPtr[K ][J ][I ] + DSonPtr[K ][J ][Ip] +
                                       DSonPtr[K ][Jp][I ] + DSonPtr[Kp][J ][I ] +
                                       DSonPtr[K ][Jp][Ip] + DSonPtr[Kp][Jp][I ] +
                                       DSonPtr[Kp][J ][Ip] + DSonPtr[Kp][Jp][Ip] );

#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
//             if father level uses fluid scheme and ELBDM_MATCH_PHASE is on, unwrap average phase to match parent phase
               if ( !amr->use_wave_flag[FaLv] && ELBDM_MATCH_PHASE ) {
                  avgphase = ELBDM_UnwrapPhase( OldPFaPtr[kk][jj][ii], avgphase );
               }

               if ( !amr->use_wave_flag[FaLv] ) {
                  if ( TVarCC & _DENS )   DFaPtr[kk][jj][ii] = avgdens;
                  if ( TVarCC & _PHAS )   PFaPtr[kk][jj][ii] = avgphase;
               } else {
#              endif
                  if ( TVarCC & _DENS )   DFaPtr[kk][jj][ii] = avgdens;
                  if ( TVarCC & _REAL )   RFaPtr[kk][jj][ii] = SQRT(avgdens) * COS(avgphase);
                  if ( TVarCC & _IMAG )   IFaPtr[kk][jj][ii] = SQRT(avgdens) * SIN(avgphase);
#              if ( ELBDM_SCHEME == ELBDM_HYBRID )
               } // if ( !amr->use_wave_flag[FaLv] ) ... else ...
#              endif
            }}}
         } // if ( ResPha )
#        endif // #if ( MODEL == ELBDM )


//       restrict the fluid data
//       ELBDM: only restrict fluid data that has not yet been restricted using phase restriction
         if ( ResFlu  &&  !Res4 ) {
         for (int v=0; v<NFluVar; v++)
         {
            const int TFluVarIdx = TFluVarIdxList[v];
            const real (*SonPtr)[PS1][PS1] = amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx];
                  real (* FaPtr)[PS1][PS1] = amr->patch[ FaFluSg][ FaLv][ FaPID]->fluid[TFluVarIdx];

            int ii, jj, kk, I, J, K, Ip, Jp, Kp;

            for (int k=0; k<PS1_half; k++)  {  K = k*2;  Kp = K+1;  kk = k + Disp_k;
            for (int j=0; j<PS1_half; j++)  {  J = j*2;  Jp = J+1;  jj = j + Disp_j;
            for (int i=0; i<PS1_half; i++)  {  I = i*2;  Ip = I+1;  ii = i + Disp_i;

               FaPtr[kk][jj][ii] = 0.125*( SonPtr[K ][J ][I ] + SonPtr[K ][J ][Ip] +
                                           SonPtr[K ][Jp][I ] + SonPtr[Kp][J ][I ] +
                                           SonPtr[K ][Jp][Ip] + SonPtr[Kp][Jp][I ] +
                                           SonPtr[Kp][J ][Ip] + SonPtr[Kp][Jp][Ip] );
            }}}
         }
         } // if ( ResFlu )


//       restrict the potential data
#        ifdef GRAVITY
         if ( ResPot )
         {
            const real (*SonPtr)[PS1][PS1] = amr->patch[SonPotSg][SonLv][SonPID]->pot;
                  real (* FaPtr)[PS1][PS1] = amr->patch[ FaPotSg][ FaLv][ FaPID]->pot;

            int ii, jj, kk, I, J, K, Ip, Jp, Kp;

            for (int k=0; k<PS1_half; k++)  {  K = k*2;  Kp = K+1;  kk = k + Disp_k;
            for (int j=0; j<PS1_half; j++)  {  J = j*2;  Jp = J+1;  jj = j + Disp_j;
            for (int i=0; i<PS1_half; i++)  {  I = i*2;  Ip = I+1;  ii = i + Disp_i;

               FaPtr[kk][jj][ii] = 0.125*( SonPtr[K ][J ][I ] + SonPtr[K ][J ][Ip] +
                                           SonPtr[K ][Jp][I ] + SonPtr[Kp][J ][I ] +
                                           SonPtr[K ][Jp][Ip] + SonPtr[Kp][Jp][I ] +
                                           SonPtr[Kp][J ][Ip] + SonPtr[Kp][Jp][Ip] );
            }}}
         } // if ( ResPot )
#        endif // #ifdef GRAVITY


//       restrict the magnetic field
//       --> currently it always works on all three B field components
//###OPTIMIZATION: coarse-grid B field on the son patch boundaries (within the same patch group)
//                 is restricted twice by two adjacent son patches
#        ifdef MHD
         if ( ResMag )
         {
            int idx_fa, idx_son0, I, J, K;

//          Bx
            const real *SonBx = amr->patch[SonMagSg][SonLv][SonPID]->magnetic[0];
                  real * FaBx = amr->patch[ FaMagSg][ FaLv][ FaPID]->magnetic[0];

            for (int k=0; k<PS1_half;   k++)  {  K = k*2;
            for (int j=0; j<PS1_half;   j++)  {  J = j*2;
            for (int i=0; i<PS1_half+1; i++)  {  I = i*2;

               const int idx_fa   = IDX321( i+Disp_i, j+Disp_j, k+Disp_k, PS1P1, PS1 );
               const int idx_son0 = IDX321( I,        J,        K,        PS1P1, PS1 );

               FaBx[idx_fa] = 0.25*( SonBx[ idx_son0                     ] +
                                     SonBx[ idx_son0 + PS1P1             ] +
                                     SonBx[ idx_son0 + PS1P1*PS1         ] +
                                     SonBx[ idx_son0 + PS1P1*PS1 + PS1P1 ] );
            }}}

//          By
            const real *SonBy = amr->patch[SonMagSg][SonLv][SonPID]->magnetic[1];
                  real * FaBy = amr->patch[ FaMagSg][ FaLv][ FaPID]->magnetic[1];

            for (int k=0; k<PS1_half;   k++)  {  K = k*2;
            for (int j=0; j<PS1_half+1; j++)  {  J = j*2;
            for (int i=0; i<PS1_half;   i++)  {  I = i*2;

               const int idx_fa   = IDX321( i+Disp_i, j+Disp_j, k+Disp_k, PS1, PS1P1 );
               const int idx_son0 = IDX321( I,        J,        K,        PS1, PS1P1 );

               FaBy[idx_fa] = 0.25*( SonBy[ idx_son0                 ] +
                                     SonBy[ idx_son0 + 1             ] +
                                     SonBy[ idx_son0 + PS1P1*PS1     ] +
                                     SonBy[ idx_son0 + PS1P1*PS1 + 1 ] );
            }}}

//          Bz
            const real *SonBz = amr->patch[SonMagSg][SonLv][SonPID]->magnetic[2];
                  real * FaBz = amr->patch[ FaMagSg][ FaLv][ FaPID]->magnetic[2];

            for (int k=0; k<PS1_half+1; k++)  {  K = k*2;
            for (int j=0; j<PS1_half;   j++)  {  J = j*2;
            for (int i=0; i<PS1_half;   i++)  {  I = i*2;

               const int idx_fa   = IDX321( i+Disp_i, j+Disp_j, k+Disp_k, PS1, PS1 );
               const int idx_son0 = IDX321( I,        J,        K,        PS1, PS1 );

               FaBz[idx_fa] = 0.25*( SonBz[ idx_son0           ] +
                                     SonBz[ idx_son0 + 1       ] +
                                     SonBz[ idx_son0 + PS1     ] +
                                     SonBz[ idx_son0 + PS1 + 1 ] );
            }}}
         } // if ( ResMag )
#        endif // ifdef MHD
      } // for (int LocalID=0; LocalID<8; LocalID++)

//    apply the same B field restriction to the data of father-sibling patches on the coarse-fine boundaries
#     ifdef MHD
      for (int s=0; s<6; s++)
      {
         const int FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#        ifdef GAMER_DEBUG
         if ( FaSibPID == -1 )
            Aux_Error( ERROR_INFO, "FaSibPID == -1 (FaLv %d, FaPID %d, s %d, SonPID0 %d) !!\n", FaLv, FaPID, s, SonPID0 );
#        endif

//       skip father patches adjacent to non-periodic boundaries
         if ( FaSibPID < -1 )    continue;

//       find the coarse-fine boundaries and copy data
         if ( amr->patch[0][FaLv][FaSibPID]->son == -1 )    MHD_CopyPatchInterfaceBField( FaLv, FaPID, s, FaMagSg );
      }
#     endif // #ifdef MHD


//    check the minimum pressure/internal energy and, when the dual-energy formalism is adopted, ensure the consistency between
//    pressure, total energy density, and the dual-energy variable
#     if ( MODEL == HYDRO  &&  !defined SRHD )
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
//       compute magnetic energy
#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( FaLv, FaPID, i, j, k, FaMagSg );
#        else
         const real Emag = NULL_REAL;
#        endif

#        ifdef DUAL_ENERGY
//       here we ALWAYS use the dual-energy variable to correct the total energy density
//       --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//           the runtime parameter DUAL_ENERGY_SWITCH here
         const bool CheckMinPres_Yes = true;
         const real UseDual2FixEngy  = HUGE_NUMBER;
         char dummy;    // we do not record the dual-energy status here

         real Passive[NCOMP_PASSIVE];
         for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = amr->patch[FaFluSg][FaLv][FaPID]->fluid[ NCOMP_FLUID + v ][k][j][i];

         Hydro_DualEnergyFix( amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMX][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMY][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMZ][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[DUAL][k][j][i],
                              Passive, dummy, EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2], CheckMinPres_Yes, MIN_PRES,
                              PassiveFloorMask, UseDual2FixEngy, Emag );

#        else // #ifdef DUAL_ENERGY

//       actually it might not be necessary to check the minimum internal energy here
         amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i]
            = Hydro_CheckMinEintInEngy( amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMX][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMY][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMZ][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i],
                                        MIN_EINT, PassiveFloorMask, Emag );
#        endif // #ifdef DUAL_ENERGY ... else ...
      } // i,j,k
#     endif // #if ( MODEL == HYDRO  &&  !defined SRHD )


//    rescale real and imaginary parts to get the correct density in ELBDM
#     if ( MODEL == ELBDM )
      if ( amr->use_wave_flag[FaLv] ) {
      real Real, Imag, Rho_Wrong, Rho_Corr, Rescale;

      if (  ( TVarCC & _DENS )  &&  ( TVarCC & _REAL )  &&  (TVarCC & _IMAG )  )
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         Real      = amr->patch[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i];
         Imag      = amr->patch[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i];
         Rho_Wrong = Real*Real + Imag*Imag;
         Rho_Corr  = amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i];

//       be careful about the negative density introduced from the round-off errors
         if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
         {
            amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i] = (real)0.0;
            Rescale = (real)0.0;
         }
         else
            Rescale = SQRT( Rho_Corr/Rho_Wrong );

         amr->patch[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i] *= Rescale;
         amr->patch[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i] *= Rescale;
      }

      } // if ( amr->use_wave_flag[FaLv] )
#     endif // # if ( MODEL == ELBDM )

   } // for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)

} // FUNCTION : Flu_FixUp_Restrict
