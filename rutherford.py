# Adam's version of ruthxx.f, original by Wu.
#
# Added Doug's "safe distance" and the distance of
# closest approach.
# Changed output format and updated fortran standards.
# Jan. 10 2008.
#
# Time values are in ns.
# VARIABLES:
# ruthi - multiplicative constant for cross section function
# The factor 1.296 converts the cross sections to millibarns (mb).
# See my physics notebook.
# ruthc - diff. crossection in c.o.m. frame
# ruthl - diff. crossection in lab frame
# theta - scat. or recoil angle, depending on loop.
# thetb - the angle for the opposite particle of "theta's" particle
# dca   - distance of closest approach in fm
# l     - multiplicative constant for dist. of closest approach
#         l is in units of fm*MeV
# dsafe - Doug's "safe" distance of closest approach
#         in fm.
#---------------------------------------------------------------------------

import matplotlib
import pylab


 92   continue


      write(*,321)
 321  format('Enter Zp,Ap (projectile), Zt,At (target) :',$)
      read *,z1,a1,z2,a2

      write(1,325)
 325  format('Output from RUTHERFORD.FOR.')
      write(1,*)''
      write(1,323)
      write(*,323)
 323  format('Zp Ap  Zt At')
      write(1,324)int(z1),int(a1),int(z2),int(a2)
      write(*,324)int(z1),int(a1),int(z2),int(a2)
 324  format(i2,1x,i3,1x,i2,1x,i3)


      cnt=49.0
      disf=12.8
      write(1,*)''
      write(1,623) 
 623  format('Set up for CHICO:')
      write(1,624) cnt
 624  format('     Angle to normal is ',f4.1,' degrees.')
      write(1,625) disf
 625  format('     Distance to detector is ',f4.1,' cm.')

        const1 = cos(cnt/radtodeg)
        const2 = sin(cnt/radtodeg)
        consta = cos((cnt+90.)/radtodeg)
        constb = sin((cnt+90.)/radtodeg)

      write(*,932)
 932  format('Enter angles--min, max, step: ',$)

      read *,amin,amax,astep

      write(*,223)amin,amax,astep

      write(1,*)''
 223  format('Angles ',f5.1,' to ',f5.1,', step ',f5.1)

      print *,' '
      write(*,922)
 922  format('Writing output to file "rutherford.dat".')
      print *,' '

      do ijkl=1,100000

        print *,'Enter beam energy, MeV (0 to exit, -1 to start over):'
        read *,e0
        IF(e0.eq.0.) CALL EXIT
        if(e0.lt.0.) goto 92
        write(*,420)e0
        write(1,420)e0
 420    format('Beam energy is ',f6.1,' MeV.')

        istep=(amax-amin)/astep+1.0

        dsafe=1.25*(a1**(1.0/3.0)+a2**(1.0/3.0)) + 5.0
        ecm=a2/(a1+a2)*e0
        am=a1*a2/(a1+a2)**2
        ruthi=(z1*z2)**2*1.296/ecm**2       !gives result in mb.
        l=1.43990*z1*z2/(2*ecm)

        thelim=180.
        if(a1.ge.a2) thelim=asin(a2/a1)*radtodeg

        write(1,423)thelim
 423    format('Maximum scattering angle in lab frame is ',
     &    f7.3,' degrees.')

        write (1,424)dsafe
 424    format('Safe distance is ',f5.2,' fm.')
        write(1,*)''


 11     format('theta',4x,'theta',4x,'theta',4x,'dist of  ',1x,
     &  'energy',1x,'dsigma/domega',1x,'dsigma/domega',1x,
     &  'dsigma/domega',1x,'beta',7x,'beta',7x,'time of flight')
 111    format('scat ',4x,'rec  ',4x,'scat ',4x,'closest  ',1x,
     &  'in lab',1x,'c.o.m.       ',1x,'laboratory   ',1x,
     &  'lab/c.o.m.   ',1x,'scat',7x,'rec ',7x,'difference    ')
 112    format('lab  ',4x,'lab  ',4x,'c.o.m.',3x,'appr (fm)',1x,
     &  '(MeV) ',1x,'(mb/sr)      ',1x,'(mb/sr)      ',1x,
     &  '(1)          ',1x,'(1) ',7x,'(1) ',7x,'(ns)          ')
        write (1,113)
 113    format('SCATTERED particle:')
        write (1,11)
        write (1,111)
        write (1,112)

          
        do 5 i=1,istep          ! Step through scattering angles in lab frame.
c         Theta is scattering angle.
          z=i-1
          theta=amin+astep*z
          if(theta.eq.0.) go to 5
          if(theta.ge.thelim) go to 45
          tcm=theta+asin(a1*sin(theta/radtodeg)/a2)*radtodeg
          thetb=(180.-tcm)/2.
          ruthc=ruthi/sin((tcm/2.)/radtodeg)**4
          e1=e0*(1.-2.*am*(1.-cos(tcm/radtodeg)))
          ratio=sin(tcm/radtodeg)**2/(cos((tcm-theta)/radtodeg)
     &      *sin(theta/radtodeg)**2)
          ruthl=ruthc*ratio
          thed1= const1*cos(theta/radtodeg) + const2*sin(theta/radtodeg)
          thed2= const1*cos(thetb/radtodeg) + const2*sin(thetb/radtodeg)
          if(theta.ge.90.)thed1= consta*cos(theta/radtodeg) 
     &      + constb*sin(theta/radtodeg)
          dis1=disf/ thed1
          dis2=disf/ thed2
          beta1=.046*sqrt(e1/a1)
          beta2=.046*sqrt((e0-e1)/a2)
          tim_s=.0333564*dis1/beta1
          tim_r=.0333564*dis2/beta2
          timdi=tim_s-tim_r
          dca=l*(1.0+(1.0/sin((tcm/2.0)/radtodeg)))


          write(1,77)theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio,beta1,
     &      beta2,timdi
 77       format(f6.2,3x,f6.2,3x,f6.2,3x,f6.1,4x,f6.1,1x,
     &          e10.3,4x,e10.3,4x,e10.3,4x,e10.3,1x,e10.3,1x,e10.3)

 5      continue
 45     if(thelim.eq.180.) go to 50

        
        write(1,117)
 117    format('SCATTERED particle, SECOND solution:')

        do 48 i=1,istep
          z=i-1
          theta=amin+astep*z
          if(theta.eq.0.) go to 48
          if(theta.ge.thelim) go to 50
          tcm=theta-asin(a1*sin(theta/radtodeg)/a2)*radtodeg
          tcm=180.-abs(tcm)
          thetb=(180.-tcm)/2.
          e1=e0*(1.-2.*am*(1.-cos(tcm/radtodeg)))
          ruthc=ruthi/sin((tcm/2.)/radtodeg)**4
          ratio=sin(tcm/radtodeg)**2/(sin(theta/radtodeg)**2
     &      *cos((tcm-theta)/radtodeg))
          ratio=abs(ratio)
          ruthl=ruthc*ratio
          dca=l*(1.0+(1.0/sin((tcm/2.0)/radtodeg)))
          write(1,77)theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio
 48     continue

 50     continue 

        write(1,126)
 126    format('RECOILING particle:')
 12     format('theta',4x,'theta',4x,'theta',4x,'dist of  ',1x,
     &  'energy',1x,'dsigma/domega',1x,'dsigma/domega',1x,
     &  'dsigma/domega')
 121    format('rec  ',4x,'scat ',4x,'rec  ',4x,'closest  ',1x,
     &  'in lab',1x,'c.o.m.       ',1x,'laboratory   ',1x,
     &  'lab/c.o.m.   ')
 122    format('lab  ',4x,'lab  ',4x,'c.o.m.',3x,'appr (fm)',1x,
     &  '(MeV) ',1x,'(mb/sr)      ',1x,'(mb/sr)      ',1x,
     &  '(1)          ')
        write (1,12)
        write (1,121)
        write (1,122)

        do 60 i=1,istep ! Step through recoil angles in lab frame.
c         Theta is now RECOIL angle.
          z=i-1
          theta=amin+astep*z
          if(theta.eq.0.) go to 60
          if(theta.ge.90.) go to 80
          e2 = e0*4.*am*cos(theta/radtodeg)**2
          tcm=acos(1.-2.*cos(theta/radtodeg)**2)*radtodeg
          thetb=atan(sin((2.*theta)/radtodeg)
     &      /(a1/a2-cos(2.*theta/radtodeg)))*radtodeg
          if(thetb.lt.0.) thetb=180.+thetb
          ruthc=ruthi/sin((tcm/2.)/radtodeg)**4
          ratio=4.*cos(theta/radtodeg)
          ruthl=ruthc*ratio
          tcm=180.-tcm
          thed1= const1*cos(theta/radtodeg) + const2*sin(theta/radtodeg)
          thed2= const1*cos(thetb/radtodeg) + const2*sin(thetb/radtodeg)
          if(thetb.ge.90.)thed1= consta*cos(thetb/radtodeg) 
     &       + constb*sin(thetb/radtodeg)
          dis1=disf/ thed1
          dis2=disf/ thed2
          beta1=.046*sqrt(e2/a2)
          beta2=.046*sqrt((e0-e2)/a1)
          tim_s=.0333564*dis1/beta1
          tim_r=.0333564*dis2/beta2
          timdi=tim_s-tim_r
          tcmscat=thetb+asin(a1*sin(thetb/radtodeg)/a2)*radtodeg
          dca=l*(1.0+(1.0/sin((tcmscat/2.0)/radtodeg)))

          write(1,77)theta,thetb,tcm,dca,e2,ruthc,ruthl,ratio

 60     continue
 80     continue

        write(1,872)
 872    format(//,128('-'),//)

        enddo  ! main loop

      end


