MODULE STAT_MODULE
    IMPLICIT NONE
    
    real       t0, t1, t2, t3, t4, t5
    save       t0, t1, t2, t3, t4, t5
!
    integer    nopx, nbx, nrorth, nitref, nrstrt
    real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv
    real       tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv
    real       tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv
    real       tmvopx, tmvbx, tgetv0, titref, trvec
    
    END MODULE STAT_MODULE