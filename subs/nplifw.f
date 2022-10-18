C*********************************************************************
C subroutine NPLot Input File Write **********************************
C            ---   -     -    -     **********************************
C Steve Gibbons Tue Dec  7 19:17:16 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out an input file for the nplot family of programs.         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LU        : Logical unit number of file.                       C
C     IWR       : Write flag. = 2 for caution, =3 for overwrite.     C
C     ISF       : =1 if a scalar file is to be included.             C
C     IRAD      : Radial node at which display defaults.             C
C     LAT       : Latitude at which display defaults.                C
C     LON       : Longitude at which display defaults.               C
C     IPLOT     : Projection flag: =1 for radial projection.         C
C                                  =2 for equatorial section.        C
C                                  =3 for meridian section.          C
C     IFN       : Function to be plotted.                            C
C                                                                    C
C                  =1 for r stream function                          C
C                  =2 for theta stream function                      C
C                  =3 for phi stream function                        C
C                                                                    C
C                  =4 for r vector component                         C
C                  =5 for theta vector component                     C
C                  =6 for phi vector component                       C
C                                                                    C
C                  =7 for r vector curl component                    C
C                  =8 for theta vector curl component                C
C                  =9 for phi vector curl component                  C
C                                                                    C
C                  =10 for helicity                                  C
C                  =11 for |v|^2  (norm of vector)                   C
C                  =12 for |w|^2  (norm of vector curl)              C
C                                                                    C
C     IAS       : Set to 1 to exclude non-axisymmetric parts.        C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NH        : Number of spherical harmonics.                     C
C     NR        : Number of equally spaced radial grid nodes.        C
C     NCL       : Number of contour levels.                          C
C     ICOL      : 0 for monochrome, =1 for colour.                   C
C     ICON      : 1 for contours to be included - 0 otherwise.       C
C     ILAB      : 1 for labels - 0 otherwise.                        C
C     ICB       : 1 for colourbar - 0 otherwise.                     C
C     IEDGE     : 1 for eq. line - 0 otherwise.                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RI        : Radial value at node 1.                            C
C     RO        : Radial value at node NR.                           C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : File name for input file.                          C
C     FNINTF    : File name for integers file.                       C
C     FNEVCF    : File name for evecs file.                          C
C     FNSCAF    : File name for scalar file.                         C
C     TITLE     : Character string for title.                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NPLIFW( LU, FNAME, IWR, FNINTF, FNEVCF, FNSCAF,
     1                   TITLE, ISF, IRAD, LAT, LON, IPLOT, IFN, IAS,
     2                   LH, NH, NR, RI, RO, NCL, ICOL, ICON, ILAB,
     3                   ICB, IEDGE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LU, IWR, ISF, IRAD, LAT, LON, IPLOT, IFN, IAS,
     1        LH, NH, NR, NCL, ICOL, ICON, ILAB, ICB, IEDGE
      CHARACTER *(*) FNAME, FNINTF, FNEVCF, FNSCAF, TITLE
      DOUBLE PRECISION RI, RO
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION ASP, LOW
      PARAMETER ( LOW = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RO.LT.LOW ) THEN
        PRINT *,' Subroutine NPLIFW.'
        PRINT *,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      ASP = RI/RO
C
      IF ( IWR.NE.2 .AND. IWR.NE.3 ) THEN
        PRINT *,' Subroutine NPLIFW.'
        PRINT *,' IWR = ', IWR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Open file
C
      CALL FOPEN ( LU, FNAME, IWR )
C
      WRITE ( LU, 781 )
      WRITE ( LU , * ) '           TITLE'
      WRITE ( LU, 783 ) TITLE
      WRITE ( LU, 782 ) 
      WRITE ( LU , * ) '           HARMONICS FILE '
      WRITE ( LU, 783 ) FNINTF
      WRITE ( LU, 782 ) 
      WRITE ( LU , * ) '           DATA FILE '
      WRITE ( LU, 783 ) FNEVCF
      IF ( ISF.EQ.1 ) THEN
        WRITE ( LU, 783 ) FNSCAF
      ELSE
        WRITE ( LU, 782 ) 
      ENDIF
      WRITE ( LU, * ) 'irad    lat     lon    iplot    ifn    ias'
      WRITE ( LU, 785 ) IRAD, LAT, LON, IPLOT, IFN, IAS
      WRITE ( LU, 782 )
      WRITE ( LU , * ) 'lh     nh     nr    inorm      ri        ro'
      WRITE ( LU , 784 ) LH, NH, NR, RI, RO, ASP
      WRITE ( LU, 782 )
      WRITE ( LU, *) 'levels      colours     contours     labels',
     1                '     colourbar     edges'
      WRITE ( LU, 785 ) NCL, ICOL, ICON, ILAB, ICB, IEDGE
      WRITE ( LU, 782 )
      WRITE ( LU, * ) ' info man_cont  cont_min  cont_max '
      WRITE ( LU, * ) '   0      0       10.0      10.0   '
      WRITE ( LU, 782 )
      WRITE ( LU, * ) 'iplot  =  projection : radial = 1; '
      WRITE ( LU, * ) ' equatorial = 2; meridional = 3. '
      WRITE ( LU, * ) 'ifn    =  fn plotted : '
      WRITE ( LU, * ) 'streamfn{r,the,phi} = 1-3;'
      WRITE ( LU, * ) 'vector compt{r,the,phi} = 4-6; '
      WRITE ( LU, * ) 'vector curl{r,the,phi} = 7-9; '
      WRITE ( LU, * ) '     helicity = 10; |v|^2 = 11; |w|^2 = 12. '
      WRITE ( LU, 781 )
C
C Close file.
C
      CALL FCLOSE ( LU, FNAME, 'Error in closing file.')
C
C Format statements
C
 781  FORMAT ('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx',
     1        'xxxxxxxxxxxxxxxxxxxx')
 782  FORMAT ('----------------------------------------------------',
     1        '--------------------')
 783  FORMAT (a)
 784  FORMAT (i4,'   ',i4,'   ',i4,'     1    ',f6.2,'     ',f6.2,
     1        '    ',f6.2)
 785  FORMAT (6(i6))
C
      RETURN
      END
C*********************************************************************
