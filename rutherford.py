#!/usr/bin/python
# Adam's version of ruthxx.f, original by Wu.
#
# Time values are in ns.
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

import math
import numpy
#import matplotlib
import matplotlib.pyplot as plt
import pylab
import sys

# Check to see if the user gave command line arguments.
print len(sys.argv), sys.argv
if len(sys.argv) == 9:
    interactive = False
    all_args = sys.argv[1:9]
    z1,a1,z2,a2,amin,amax,astep,e0 = [float(x) for x in all_args]
else:
    print 'Command-line usage:\nrutherford.py Zp Ap Zt At theta_min theta_max theta_step beam_energy'
    print '  e.g. rutherford.py 54 136 74 186 10 170 1 800 > out.txt'
    print '                     Zp  Ap Zt  At |angles| Ebeam [MeV] (redirected to file "out.txt".'
    interactive = True
# cnt is the angle of a normal from the target position to a planar detector.
# disf is the distance along that normal from target to detector.

cnt  = 49.0  # degrees
disf = 12.8  # cm

# Constants for the time-of-flight calculation to a plane detector.
const1 = math.cos(math.radians(cnt))
const2 = math.sin(math.radians(cnt))
consta = math.cos(math.radians(cnt + 90.))
constb = math.sin(math.radians(cnt + 90.))

def get_floats(prompt_string = ''):
    temp_string = raw_input(prompt_string)
    temp_string = temp_string.replace(',',' ')
    temp_list   = temp_string.split()
    return [float(i) for i in temp_list]

if interactive:
    z1,a1,z2,a2 = get_floats('Enter Zp,Ap (projectile), Zt,At (target) : ')
    amin,amax,astep = get_floats('Enter angles--min, max, step: ')
    e0 = float(raw_input('Enter beam energy in MeV: '))

# Reprint parameters.
for one in ['z1','a1','z2','a2','amin','amax','astep','e0']:
    print one.ljust(6) + ": " + str(eval(one))

print 'Time-of-flight difference will be calculated for CHICO:\nAngle to normal is ' + str(cnt) + ' degrees.\nDistance to detector is ' + str(disf) + ' degrees.'

# Lists for d-TOF plot.
theta_list = []
timdi_list = []
recoil_theta_list = []
recoil_timdi_list = []

dsafe=1.25*(a1**(1.0/3.0)+a2**(1.0/3.0)) + 5.0  # Doug's safe criterion.
ecm=a2/(a1+a2)*e0
am=a1*a2/(a1+a2)**2
ruthi=(z1*z2)**2*1.296/ecm**2       # gives result in mb.
l=1.43990*z1*z2/(2*ecm)

thelim=180.
if a1 >= a2:
    thelim = math.radians(math.asin(a2/a1))

print 'Maximum scattering angle in lab frame is ' + str(thelim) + ' degrees.'

print 'Safe distance is ' + str(dsafe) + ' fm.'

header = [
        'theta'.rjust(6) + 'theta'.rjust(8) + 'theta'.rjust(8)  + 'dist of'.rjust(11)   +  'energy'.rjust(8) + 'dsigma/domega'.rjust(15) + 'dsigma/domega'.rjust(15) +  'dsigma/domega'.rjust(15) + 'beta'.rjust(8) + 'beta'.rjust(8) + 'time of flight'.rjust(16),
        'scat '.rjust(6) + 'rec  '.rjust(8) + 'scat'.rjust(8)   + 'closest'.rjust(11)   +  'in lab'.rjust(8) + 'c.o.m.'.rjust(15)        + 'laboratory'.rjust(15)    +  'lab/c.o.m.'.rjust(15)    + 'scat'.rjust(8) + 'rec'.rjust(8)  + 'difference'.rjust(16),
        'lab  '.rjust(6) + 'lab  '.rjust(8) + 'c.o.m.'.rjust(8) + 'appr (fm)'.rjust(11) +  '(MeV)'.rjust(8)  + '(mb/sr)'.rjust(15)       + '(mb/sr)'.rjust(15)       +  '(1)'.rjust(15)           + '(1)'.rjust(8)  + '(1)'.rjust(8)  + '(ns)'.rjust(16)
         ]

print 'SCATTERED particle:'
for line in header:
    print line

for theta in numpy.linspace(amin, amax, int((amax - amin) / astep) + 1): # do 5 i=1,istep          ! Step through scattering angles in lab frame.
    # Theta is scattering angle.
    if theta >= thelim:
        break

    tcm = theta + math.degrees(math.asin(a1*math.sin(math.radians(theta))/a2))
    thetb = (180.-tcm)/2.
    ruthc = ruthi/math.sin(math.radians(tcm/2.))**4
    e1 = e0*(1.-2.*am*(1.-math.cos(math.radians(tcm))))
    ratio = math.sin(math.radians(tcm))**2/(math.cos(math.radians(tcm-theta)) *math.sin(math.radians(theta))**2)
    ruthl = ruthc*ratio
    thed1 =  const1*math.cos(math.radians(theta)) + const2*math.sin(math.radians(theta))
    thed2 =  const1*math.cos(math.radians(thetb)) + const2*math.sin(math.radians(thetb))
    if theta >= 90.0:
        thed1 =  consta*math.cos(math.radians(theta)) + constb*math.sin(math.radians(theta))
    dis1 = disf/ thed1
    dis2 = disf/ thed2
    beta1 = .046*math.sqrt(e1/a1)
    beta2 = .046*math.sqrt((e0-e1)/a2)
    tim_s = .0333564*dis1/beta1
    tim_r = .0333564*dis2/beta2
    timdi = tim_s-tim_r

    # Add to lists for delta-TOF plot:
    theta_list.append(theta)
    timdi_list.append(timdi)

    dca = l*(1.0+(1.0/math.sin(math.radians(tcm/2.0))))

    print '{0:6.2f}{1:8.2f}{2:8.2f}{3:11.2f}{4:8.2f}{5:15.2e}{6:15.2e}{7:15.2e}{8:8.3f}{9:8.3f}{10:16.2f}'.format(theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio,beta1,beta2,timdi)


if thelim < 179.0:

    header = [
            'theta'.rjust(6) + 'theta'.rjust(8) + 'theta'.rjust(8)  + 'dist of'.rjust(11)   +  'energy'.rjust(8) + 'dsigma/domega'.rjust(15) + 'dsigma/domega'.rjust(15) +  'dsigma/domega'.rjust(15) + 'beta'.rjust(8) + 'beta'.rjust(8) + 'time of flight'.rjust(16),
            'scat '.rjust(6) + 'rec  '.rjust(8) + 'scat'.rjust(8)   + 'closest'.rjust(11)   +  'in lab'.rjust(8) + 'c.o.m.'.rjust(15)        + 'laboratory'.rjust(15)    +  'lab/c.o.m.'.rjust(15)    + 'scat'.rjust(8) + 'rec'.rjust(8)  + 'difference'.rjust(16),
            'lab  '.rjust(6) + 'lab  '.rjust(8) + 'c.o.m.'.rjust(8) + 'appr (fm)'.rjust(11) +  '(MeV)'.rjust(8)  + '(mb/sr)'.rjust(15)       + '(mb/sr)'.rjust(15)       +  '(1)'.rjust(15)           + '(1)'.rjust(8)  + '(1)'.rjust(8)  + '(ns)'.rjust(16)
             ]

    print 'SCATTERED particle, SECOND solution:'

    for line in header:
        print line

    for theta in numpy.linspace(amin, thelim, int((thelim - amin) / astep) + 1):

        if theta == 0.0:
            continue
        if theta >= thelim:
            break
        tcm=theta-math.degrees(math.asin(a1*math.sin(math.radians(theta))/a2))
        tcm=180.-abs(tcm)
        thetb=(180.-tcm)/2.
        e1=e0*(1.-2.*am*(1.-math.cos(math.radians(tcm))))
        ruthc=ruthi/math.sin(math.radians(tcm/2.))**4
        ratio=math.sin(math.radians(tcm))**2/(math.sin(math.radians(theta))**2
        *math.cos(math.radians(tcm-theta)))
        ratio=abs(ratio)
        ruthl=ruthc*ratio
        dca=l*(1.0+(1.0/math.sin(math.radians(tcm/2.0))))

        print '{0:6.2f}{1:8.2f}{2:8.2f}{3:11.2f}{4:8.2f}{5:15.2e}{6:15.2e}{7:15.2e}{8:8.3f}{9:8.3f}{10:16.2f}'.format(theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio)

#     12 format('theta',4x,'theta',4x,'theta',4x,'dist of  ',1x,'energy',1x,'dsigma/domega',1x,'dsigma/domega',1x,'dsigma/domega')
#     121format('rec  ',4x,'scat ',4x,'rec  ',4x,'closest  ',1x,'in lab',1x,'c.o.m.       ',1x,'laboratory   ',1x,'lab/c.o.m.   ')
#     122format('lab  ',4x,'lab  ',4x,'c.o.m.',3x,'appr (fm)',1x,'(MeV) ',1x,'(mb/sr)      ',1x,'(mb/sr)      ',1x,'(1)          ')

# Recoiling target.
header = [
        'theta'.rjust(6) + 'theta'.rjust(8) + 'theta'.rjust(8)  + 'dist of'.rjust(11)   +  'energy'.rjust(8) + 'dsigma/domega'.rjust(15) + 'dsigma/domega'.rjust(15) +  'dsigma/domega'.rjust(15) + 'beta'.rjust(8) + 'beta'.rjust(8) + 'time of flight'.rjust(16),
        'rec'.rjust(6) + 'scat'.rjust(8) + 'rec'.rjust(8)   + 'closest'.rjust(11)   +  'in lab'.rjust(8) + 'c.o.m.'.rjust(15)        + 'laboratory'.rjust(15)    +  'lab/c.o.m.'.rjust(15)    + 'scat'.rjust(8) + 'rec'.rjust(8)  + 'difference'.rjust(16),
        'lab  '.rjust(6) + 'lab  '.rjust(8) + 'c.o.m.'.rjust(8) + 'appr (fm)'.rjust(11) +  '(MeV)'.rjust(8)  + '(mb/sr)'.rjust(15)       + '(mb/sr)'.rjust(15)       +  '(1)'.rjust(15)           + '(1)'.rjust(8)  + '(1)'.rjust(8)  + '(ns)'.rjust(16)
         ]

print 'RECOILING particle:'

for line in header:
    print line

for theta in numpy.linspace(amin, min(amax,89.), int((min(amax,89.) - amin) / astep) + 1):
    # theta is now RECOIL angle.
    e2 = e0*4.*am*math.cos(math.radians(theta))**2
    tcm=math.degrees(math.acos(1.-2.*math.cos(math.radians(theta))**2))
    thetb=math.degrees(math.atan(math.sin(math.radians(2.*theta)) /(a1/a2-math.cos(2.*math.radians(theta)))))
    if thetb < 0.0:
        thetb=180.+thetb
    ruthc=ruthi/math.sin(math.radians(tcm/2.))**4
    ratio=4.*math.cos(math.radians(theta))
    ruthl=ruthc*ratio
    tcm=180.-tcm
    thed1= const1*math.cos(math.radians(theta)) + const2*math.sin(math.radians(theta))
    thed2= const1*math.cos(math.radians(thetb)) + const2*math.sin(math.radians(thetb))
    if(thetb >= 90.0):
        thed1= consta*math.cos(math.radians(thetb)) + constb*math.sin(math.radians(thetb))
    dis1=disf/ thed1
    dis2=disf/ thed2
    beta1=.046*math.sqrt(e2/a2)
    beta2=.046*math.sqrt((e0-e2)/a1)
    tim_s=.0333564*dis1/beta1
    tim_r=.0333564*dis2/beta2
    timdi=tim_s-tim_r
    tcmscat=thetb+math.degrees(math.asin(a1*math.sin(math.radians(thetb))/a2))
    dca=l*(1.0+(1.0/math.sin(math.radians(tcmscat/2.0))))

    print '{0:6.2f}{1:8.2f}{2:8.2f}{3:11.2f}{4:8.2f}{5:15.2e}{6:15.2e}{7:15.2e}{8:8.3f}{9:8.3f}{10:16.2f}'.format(theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio,beta1,beta2,timdi)

    # Add to lists for delta-TOF plot:
    recoil_theta_list.append(theta)
    recoil_timdi_list.append(timdi)


proj_plot = plt.plot(theta_list, timdi_list)
rec_plot  = plt.plot(recoil_theta_list, recoil_timdi_list)

# Use "proxy artist" (dummy, not plotted) to add a legend.
proxy_p = plt.Rectangle((0, 0), 1, 1, fc="b")
proxy_r = plt.Rectangle((0, 0), 1, 1, fc="g")
plt.legend([proxy_p, proxy_r], ["projectile","target"])

plt.xlabel('scatter / recoil angle (deg)')
plt.ylabel('time-of-flight difference (ns)')
pylab.show()

