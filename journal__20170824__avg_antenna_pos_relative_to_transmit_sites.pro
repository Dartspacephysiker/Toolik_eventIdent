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
PRO JOURNAL__20170824__AVG_ANTENNA_POS_RELATIVE_TO_TRANSMIT_SITES

  COMPILE_OPT IDL2,STRICTARRSUBS

  rel_to_centroid  = 1

  ftToM            = 0.3048D
  ftTokm           = 0.0003048D

  antsForAvg       = ["West","South","Center"]
  ;; ants = {name : ["West", $
  ;;                 "South", $
  ;;                 "Center"], $
  ;;         lat  : [68.626134D, $
  ;;                 68.625536D, $
  ;;                 68.625718D], $
  ;;         lon  : [-149.593213D, $
  ;;                -149.593789D, $
  ;;                -149.592700D], $
  ;;         alt  : [2513.7D, $
  ;;                 2506.0D, $
  ;;                 2510.0D] * ftTokm}
  ants = {name : ["Polarization", $
                  "West", $
                  "South", $
                  "Center", $
                  "Near"], $
          lat  : [68.625956D, $
                  68.626134D, $
                  68.625536D, $
                  68.625718D, $
                  68.626300D], $
          lon  : [-149.594232D, $
                  -149.593213D, $
                  -149.593789D, $
                  -149.592700D, $
                  -149.594680D], $
          alt  : [2516.3D, $
                  2513.7D, $
                  2506.0D, $
                  2510.0D, $
                  2482.5D] * ftTokm}

  sites = {name : ["HillWest", $
                  "SlopeSW"], $
           lat  : [68.D + 37.517D/60.D, $
                   68.D + 37.177D/60.D], $
           lon  : [-149.D - 36.839D/60.D, $
                   -149.D - 36.414D/60.D], $
           alt  : [2600.D, $
                   2600.D] * ftTokm}

  nAntenna         = N_ELEMENTS(ants.name)
  nSite            = N_ELEMENTS(sites.name)

  tmplt = {name  : '', $
           lat   : 0.D, $
           lon   : 0.D, $
           alt   : 0.D, $
           ECEF  : {abs : MAKE_ARRAY(3,VALUE=0.D,/DOUBLE), $
                    rel : MAKE_ARRAY(3,VALUE=0.D,/DOUBLE), $
                    rel_dist : 0.D}, $
           local : {abs : MAKE_ARRAY(3,VALUE=0.D,/DOUBLE), $
                    rel : MAKE_ARRAY(3,VALUE=0.D,/DOUBLE), $
                    rel_dist : 0.D}, $
           rel_name : ''}

  wrt_datOne_i   = (WHERE(STRUPCASE(ants.name) EQ 'SOUTH'))[0]
  tmplt.rel_name = ants.name[wrt_datOne_i]

  ant   = !NULL
  FOR k=0,nAntenna-1 DO BEGIN
     tmpStr      = tmplt
     tmpStr.name = ants.name[k]
     tmpStr.lat  = ants.lat[k]
     tmpStr.lon  = ants.lon[k]
     tmpStr.alt  = ants.alt[k]
     ant         = [ant, TEMPORARY(tmpStr)]
  ENDFOR

  site  = !NULL
  FOR k=0,nSite-1 DO BEGIN
     tmpStr      = tmplt
     tmpStr.name = sites.name[k]
     tmpStr.lat  = sites.lat[k]
     tmpStr.lon  = sites.lon[k]
     tmpStr.alt  = sites.alt[k]
     site        = [site, TEMPORARY(tmpStr)]
  ENDFOR

  ;; ECEF             = MAKE_ARRAY(nAntenna,3,VALUE=0.0D,/DOUBLE)
  ;; local            = ECEF

  ;; relatif_ECEF     = ECEF
  ;; relatif_local    = ECEF

  ;; mags_ECEF        = ECEF
  ;; mags_local       = ECEF

  ;;Get unit vecs
  tmpLambda               = ant[wrt_datOne_i].lon*!DTOR
  tmpPhi                  = ant[wrt_datOne_i].lat*!DTOR

  ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
  GEOPACK_GEODGEO_08,ant[wrt_datOne_i].alt,tmpPhi, $
                     refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
  ;; PRINT,refAlt_GEO

  GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda, $
                    refAlt_x,refAlt_y,refAlt_z,/TO_RECT

  ;; IF KEYWORD_SET(wrt_datOne_i) THEN BEGIN
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

  ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
  GEOPACK_GEODGEO_08,ant[wrt_datOne_i].alt,tmpPhi, $
                     refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
  ;; PRINT,refAlt_GEO

  GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda, $
                    refAlt_x,refAlt_y,refAlt_z,/TO_RECT

  PRINT,ant[wrt_datOne_i].name,refAlt_x,refAlt_y,refAlt_z

  rel_ECEF  = [refAlt_x,refAlt_y,refAlt_z]*1000.D
  rel_local = [DOTP(REFORM(rel_ECEF),eEast_Car_GEO), $
               DOTP(REFORM(rel_ECEF),eNorth_Car_GEO), $ 
               DOTP(REFORM(rel_ECEF),eUp_Car_GEO)]
  
  FOR i=0,nAntenna-1 DO BEGIN

     tmpLambda               = ant[i].lon*!DTOR
     tmpPhi                  = ant[i].lat*!DTOR

     ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
     GEOPACK_GEODGEO_08,ant[i].alt,tmpPhi, $
                        refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
     ;; PRINT,refAlt_GEO

     GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda, $
                       refAlt_x,refAlt_y,refAlt_z,/TO_RECT

     ;; PRINT,ant[i].name,refAlt_x,refAlt_y,refAlt_z

     ant[i].ECEF.abs        = [refAlt_x,refAlt_y,refAlt_z]*1000.D

     ant[i].local.abs[0]    = DOTP(REFORM(ant[i].ECEF.abs),eEast_Car_GEO)
     ant[i].local.abs[1]    = DOTP(REFORM(ant[i].ECEF.abs),eNorth_Car_GEO)
     ant[i].local.abs[2]    = DOTP(REFORM(ant[i].ECEF.abs),eUp_Car_GEO)

     ant[i].ECEF.rel        = ant[i].ECEF.abs - rel_ECEF
     ant[i].local.rel       = ant[i].local.abs - rel_local

     ant[i].ECEF.rel_dist   = SQRT(ant[i].ECEF.rel[0]*ant[i].ECEF.rel[0] + $
                                   ant[i].ECEF.rel[1]*ant[i].ECEF.rel[1] + $
                                   ant[i].ECEF.rel[2]*ant[i].ECEF.rel[2])
     ant[i].local.rel_dist  = SQRT(ant[i].local.rel[0]*ant[i].local.rel[0] + $
                                   ant[i].local.rel[1]*ant[i].local.rel[1] + $
                                   ant[i].local.rel[2]*ant[i].local.rel[2])

     ;; PRINT,antennas[i],",", $
     ;;       ant[i].ECEF.rel[0],",", $
     ;;       ant[i].ECEF.rel[1],",", $
     ;;       ant[i].ECEF.rel[2], $
     ;;       ant[i].ECEF.rel_dist

  ENDFOR

  ;;calc avg stuff
  antsforAvg_i           = !NULL
  FOREACH avgAnt,antsForAvg DO BEGIN
     antsForAvg_i        = [antsForAvg_i,(WHERE(STRUPCASE(ants.name) EQ STRUPCASE(avgAnt)))[0]]
  ENDFOREACH
  
  avgStr                 = tmplt
  avgStr.name            = 'Avg (So-Ce-W)'
  avgStr.ECEF.abs        = MEAN(ant[antsForAvg_i].ECEF.abs,DIMENSION=2)
  avgStr.local.abs       = MEAN(ant[antsForAvg_i].local.abs,DIMENSION=2)
  avgStr.ECEF.rel        = avgStr.ECEF.abs - rel_ECEF
  avgStr.local.rel       = avgStr.local.abs - rel_local

  avgStr.ECEF.rel_dist   = SQRT(avgStr.ECEF.rel[0]*avgStr.ECEF.rel[0] + $
                                avgStr.ECEF.rel[1]*avgStr.ECEF.rel[1] + $
                                avgStr.ECEF.rel[2]*avgStr.ECEF.rel[2])
  avgStr.local.rel_dist  = SQRT(avgStr.local.rel[0]*avgStr.local.rel[0] + $
                                avgStr.local.rel[1]*avgStr.local.rel[1] + $
                                avgStr.local.rel[2]*avgStr.local.rel[2])

  GEOPACK_SPHCAR_08,avgStr.ECEF.abs[0]/1000.D,avgStr.ECEF.abs[1]/1000.D,avgStr.ECEF.abs[2]/1000.D, $
                    refAlt_GEO,refTheta_GEO,tmpLambda,/TO_SPH

  ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
  GEOPACK_GEODGEO_08,refAlt_GEO,refTheta_GEO, $
                     altitude,tmpPhi,/TO_GEODETIC

  avgStr.lon             = tmpLambda * 180.D / !PI - 360.D
  avgStr.lat             = tmpPhi * 180.D / !PI
  avgStr.alt             = altitude

  ant                    = [ant,avgStr]
  nAntenna++

  ;;Make the average the relative loc
  IF KEYWORD_SET(rel_to_centroid) THEN BEGIN

     rel_ECEF   = ant[-1].ECEF.abs
     rel_local  = ant[-1].local.abs

     ;;Now redo locs
     FOR i=0,nAntenna-1 DO BEGIN

        tmpLambda               = ant[i].lon*!DTOR
        tmpPhi                  = ant[i].lat*!DTOR

        ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
        GEOPACK_GEODGEO_08,ant[i].alt,tmpPhi, $
                           refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
        ;; PRINT,refAlt_GEO

        GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda, $
                          refAlt_x,refAlt_y,refAlt_z,/TO_RECT

        ;; PRINT,ant[i].name,refAlt_x,refAlt_y,refAlt_z

        ant[i].ECEF.abs        = [refAlt_x,refAlt_y,refAlt_z]*1000.D

        ant[i].local.abs[0]    = DOTP(REFORM(ant[i].ECEF.abs),eEast_Car_GEO)
        ant[i].local.abs[1]    = DOTP(REFORM(ant[i].ECEF.abs),eNorth_Car_GEO)
        ant[i].local.abs[2]    = DOTP(REFORM(ant[i].ECEF.abs),eUp_Car_GEO)

        ant[i].ECEF.rel        = ant[i].ECEF.abs - rel_ECEF
        ant[i].local.rel       = ant[i].local.abs - rel_local

        ant[i].ECEF.rel_dist   = SQRT(ant[i].ECEF.rel[0]*ant[i].ECEF.rel[0] + $
                                      ant[i].ECEF.rel[1]*ant[i].ECEF.rel[1] + $
                                      ant[i].ECEF.rel[2]*ant[i].ECEF.rel[2])
        ant[i].local.rel_dist  = SQRT(ant[i].local.rel[0]*ant[i].local.rel[0] + $
                                      ant[i].local.rel[1]*ant[i].local.rel[1] + $
                                      ant[i].local.rel[2]*ant[i].local.rel[2])

        ;; PRINT,antennas[i],",", $
        ;;       ant[i].ECEF.rel[0],",", $
        ;;       ant[i].ECEF.rel[1],",", $
        ;;       ant[i].ECEF.rel[2], $
        ;;       ant[i].ECEF.rel_dist

     ENDFOR

  ENDIF

  ;;Now sites
  FOR i=0,nSite-1 DO BEGIN

     tmpLambda               = site[i].lon*!DTOR
     tmpPhi                  = site[i].lat*!DTOR

     ;;**Get colatitude (refTheta) and altitude in GEO coordinates for geodetic reference altitude (sea level)
     GEOPACK_GEODGEO_08,site[i].alt,tmpPhi, $
                        refAlt_GEO,refTheta_GEO,/TO_GEOCENTRIC ;/TO_GEODETIC (see GEOPACK_2008 documentation, or type GEOPACK_HELP)
     ;; PRINT,refAlt_GEO

     GEOPACK_SPHCAR_08,refAlt_GEO,refTheta_GEO,tmpLambda, $
                       refAlt_x,refAlt_y,refAlt_z,/TO_RECT

     ;; PRINT,site[i].name,refAlt_x,refAlt_y,refAlt_z

     site[i].ECEF.abs        = [refAlt_x,refAlt_y,refAlt_z]*1000.D

     site[i].local.abs[0]    = DOTP(REFORM(site[i].ECEF.abs),eEast_Car_GEO)
     site[i].local.abs[1]    = DOTP(REFORM(site[i].ECEF.abs),eNorth_Car_GEO)
     site[i].local.abs[2]    = DOTP(REFORM(site[i].ECEF.abs),eUp_Car_GEO)

     site[i].ECEF.rel        = site[i].ECEF.abs - rel_ECEF
     site[i].local.rel       = site[i].local.abs - rel_local

     site[i].ECEF.rel_dist   = SQRT(site[i].ECEF.rel[0]*site[i].ECEF.rel[0] + $
                                   site[i].ECEF.rel[1]*site[i].ECEF.rel[1] + $
                                   site[i].ECEF.rel[2]*site[i].ECEF.rel[2])
     site[i].local.rel_dist  = SQRT(site[i].local.rel[0]*site[i].local.rel[0] + $
                                   site[i].local.rel[1]*site[i].local.rel[1] + $
                                   site[i].local.rel[2]*site[i].local.rel[2])

  ENDFOR


  PRINT,""
  PRINT,"****"
  PRINT,"LOCAL"
  PRINT,"****"
  PRINT,FORMAT='(A0,T15,A0,T25,A0,T35,A0,T45,A0,T60,A0)', $
        'Antenna',"East (m)","North (m)","Up (m)","Bearing (deg)","Dist (m)"
  FOR i=0,nAntenna-1 DO BEGIN

     PRINT,FORMAT='(A0,T15,F0.2,T25,F0.2,T35,F0.2,T45,F0.2,T60,F0.2)', $
           ant[i].name, $
           ant[i].local.rel[0], $
           ant[i].local.rel[1], $
           ant[i].local.rel[2], $
           ATAN(ant[i].local.rel[1],ant[i].local.rel[0])*!RADEG, $
           ant[i].local.rel_dist

  ENDFOR

  PRINT,""
  PRINT,FORMAT='(A0,T15,A0,T25,A0,T35,A0,T45,A0,T60,A0)', $
        'Site',"East (m)","North (m)","Up (m)","Bearing (deg)","Dist (m)"
  FOR i=0,nSite-1 DO BEGIN

     PRINT,FORMAT='(A0,T15,F0.2,T25,F0.2,T35,F0.2,T45,F0.2,T60,F0.2)', $
           site[i].name, $
           site[i].local.rel[0], $
           site[i].local.rel[1], $
           site[i].local.rel[2], $
           ATAN(site[i].local.rel[1],site[i].local.rel[0])*!RADEG, $
           site[i].local.rel_dist

  ENDFOR

;;IDL> -176.14+90
;;      -86.139999
;;IDL> -176.14-90
;;      -266.14001
;;IDL> (-176.14-90) MOD 180
;;      -86.140015
;;IDL> ((-176.14-90) MOD 180) + 90
;;       3.8599854
;;IDL> ((-176.14-90) MOD 180) + 180
;;       93.859985
;;IDL> ((-128.91-90) MOD 180) + 180
;;       141.09000
;;IDL> 


  ;;Mags check out?
  ;; PRINT,DOTP(eNorth_Car_GEO,eEast_Car_GEO)
  ;; PRINT,DOTP(eNorth_Car_GEO,eUp_Car_GEO)
  ;; PRINT,DOTP(eEast_Car_GEO,eUp_Car_GEO)
  ;; PRINT,VECNORM(eNorth_Car_GEO)
  ;; PRINT,VECNORM(eEast_Car_GEO)
  ;; PRINT,VECNORM(eUp_Car_GEO)

  PRINT,''
  PRINT,FORMAT='(A0,T15,A0,T30,A0,T50,A0)', $
        'Antenna',"Lat (deg)","Long (deg)","Alt (m)"
  FOR i=0,nAntenna-1 DO BEGIN

     PRINT,FORMAT='(A0,T15,F0.5,T30,F0.5,T50,F0.5)', $
           ant[i].name, $
           ant[i].lat, $
           ant[i].lon, $
           ant[i].alt*1000.D

  ENDFOR

  PRINT,''
  PRINT,FORMAT='(A0,T15,A0,T30,A0,T50,A0)', $
        'Site',"Lat (deg)","Long (deg)","Alt (m)"
  FOR i=0,nSite-1 DO BEGIN

     PRINT,FORMAT='(A0,T15,F0.5,T30,F0.5,T50,F0.5)', $
           site[i].name, $
           site[i].lat, $
           site[i].lon, $
           site[i].alt*1000.

  ENDFOR


END

