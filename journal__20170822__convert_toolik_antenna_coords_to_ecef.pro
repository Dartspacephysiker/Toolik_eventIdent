;2017/08/23
FUNCTION DOTP,v1,v2
  RETURN,(TRANSPOSE(v1) # v2)[0]
END
FUNCTION VECNORM,vec
  RETURN,(SQRT(TRANSPOSE(vec) # vec))[0]
END
FUNCTION VNORMALIZE,vec
  ;; RETURN,[vec[0],vec[1],vec[2]]/SQRT(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])
  RETURN,([vec[0],vec[1],vec[2]]/VECNORM(vec))
END
FUNCTION CROSSP_NORMED,v1,v2
  tmp = CROSSP(v1,v2)
  RETURN,VNORMALIZE(tmp)
END
PRO JOURNAL__20170822__CONVERT_TOOLIK_ANTENNA_COORDS_TO_ECEF

  COMPILE_OPT IDL2,STRICTARRSUBS

  ftToM            = 0.3048D
  ftTokm           = 0.0003048D

  antennas         = ["Polarization", $
                      "West", $
                      "South", $
                      "Center", $
                      "Near"]
  
  lats             = [68.625956D, $
                      68.626134D, $
                      68.625536D, $
                      68.625718D, $
                      68.626300D]
  
  lons             = [-149.594232D, $
                      -149.593213D, $
                      -149.593789D, $
                      -149.592700D, $
                      -149.594680D]
  
  alts             = [2516.3D, $
                      2513.7D, $
                      2506.0D, $
                      2510.0D, $
                      2482.5D] * ftTokm
  ;; alts             = REPLICATE(MEAN(alts),N_ELEMENTS(alts))

  nAntenna         = N_ELEMENTS(antennas)

  ECEF             = MAKE_ARRAY(nAntenna,3,VALUE=0.0D,/DOUBLE)
  local            = ECEF

  relatif_ECEF     = ECEF
  relatif_local    = ECEF

  mags_ECEF        = ECEF
  mags_local       = ECEF
  wrt_datOne_i     = (WHERE(STRUPCASE(antennas) EQ 'SOUTH'))[0]

  PRINT,"Antenna,ECEF_x,ECEF_y,ECEF_z"
  FOR i=0,nAntenna-1 DO BEGIN

     tmpLambda               = lons[i]*!DTOR
     tmpPhi                  = lats[i]*!DTOR

     ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
     GEOPACK_GEODGEO_08,alts[i],tmpPhi, $
                        refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
     ;; PRINT,refAlt_GEO

     GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda, $
                       refAlt_x,refAlt_y,refAlt_z,/TO_RECT

     IF i EQ wrt_datOne_i THEN BEGIN
        GEOPACK_BSPCAR_08,refTheta_GEO,tmpLambda, $
                          0.D,0.D,1.D, $
                          eEast_x_GEO,eEast_y_GEO,eEast_z_GEO
        eEast_Car_GEO    = [TEMPORARY(eEast_x_GEO),TEMPORARY(eEast_y_GEO),TEMPORARY(eEast_z_GEO)]

        GEOPACK_BSPCAR_08,refTheta_GEO,tmpLambda, $
                          0.D,-1.D,0.D, $
                          eNorth_x_GEO,eNorth_y_GEO,eNorth_z_GEO
        eNorth_Car_GEO   = [TEMPORARY(eNorth_x_GEO),TEMPORARY(eNorth_y_GEO),TEMPORARY(eNorth_z_GEO)]

        GEOPACK_BSPCAR_08,refTheta_GEO,tmpLambda, $
                          1.D,0.D,0.D, $
                          eUp_x_GEO,eUp_y_GEO,eUp_z_GEO
        eUp_Car_GEO   = [TEMPORARY(eUp_x_GEO),TEMPORARY(eUp_y_GEO),TEMPORARY(eUp_z_GEO)]

     ENDIF

     PRINT,antennas[i],refAlt_x,refAlt_y,refAlt_z

     ECEF[i,*] = [refAlt_x,refAlt_y,refAlt_z]

  ENDFOR

  ;;Mags check out?
  PRINT,DOTP(eNorth_Car_GEO,eEast_Car_GEO)
  PRINT,DOTP(eNorth_Car_GEO,eUp_Car_GEO)
  PRINT,DOTP(eEast_Car_GEO,eUp_Car_GEO)
  PRINT,VECNORM(eNorth_Car_GEO)
  PRINT,VECNORM(eEast_Car_GEO)
  PRINT,VECNORM(eUp_Car_GEO)

  relatif_ECEF = ECEF
  FOR k=0,2 DO relatif_ECEF[*,k] -= REPLICATE(ECEF[wrt_datOne_i,k],nAntenna)

  relatif_ECEF *= 1000.D

  PRINT,''
  PRINT,"****"
  PRINT,"ECEF"
  PRINT,"****"
  FOR i=0,nAntenna-1 DO BEGIN

        local[i,0] = DOTP(REFORM(ECEF[i,*]),eNorth_Car_GEO)
        local[i,1] = DOTP(REFORM(ECEF[i,*]),eEast_Car_GEO)
        local[i,2] = DOTP(REFORM(ECEF[i,*]),eUp_Car_GEO)

     mags_ECEF[i] = SQRT(relatif_ECEF[i,0]*relatif_ECEF[i,0] + $
                         relatif_ECEF[i,1]*relatif_ECEF[i,1] + $
                         relatif_ECEF[i,2]*relatif_ECEF[i,2])

     PRINT,antennas[i],",", $
           relatif_ECEF[i,0],",", $
           relatif_ECEF[i,1],",", $
           relatif_ECEF[i,2], $
           mags_ECEF[i]

  ENDFOR

  relatif_local = local
  FOR k=0,2 DO relatif_local[*,k] -= REPLICATE(local[wrt_datOne_i,k],nAntenna)

  relatif_local *= 1000.D

  PRINT,""
  PRINT,"****"
  PRINT,"LOCAL"
  PRINT,"****"
  PRINT,FORMAT='(A0,T15,A0,T25,A0,T35,A0,T45,A0,T60,A0)', $
        'Antenna',"X (m)","Y (m)","Z (m)","Bearing (deg)","Dist (m)"
  FOR i=0,nAntenna-1 DO BEGIN

        mags_local[i] = SQRT(relatif_local[i,0]*relatif_local[i,0] + $
                             relatif_local[i,1]*relatif_local[i,1] + $
                             relatif_local[i,2]*relatif_local[i,2])

        PRINT,FORMAT='(A0,T15,F0.2,T25,F0.2,T35,F0.2,T45,F0.2,T60,F0.2)', $
              antennas[i], $
              relatif_local[i,0], $
              relatif_local[i,1], $
              relatif_local[i,2], $
              ATAN(relatif_local[i,0],relatif_local[i,1])*!RADEG, $
              mags_local[i]

  ENDFOR

END
