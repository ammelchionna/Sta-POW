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

# 11 format('theta',4x,'theta',4x,'theta',4x,'dist of  ',1x,
#    'energy',1x,'dsigma/domega',1x,'dsigma/domega',1x,
#    'dsigma/domega',1x,'beta',7x,'beta',7x,'time of flight')
# 111format('scat ',4x,'rec  ',4x,'scat ',4x,'closest  ',1x,
#    'in lab',1x,'c.o.m.       ',1x,'laboratory   ',1x,
#    'lab/c.o.m.   ',1x,'scat',7x,'rec ',7x,'difference    ')
# 112format('lab  ',4x,'lab  ',4x,'c.o.m.',3x,'appr (fm)',1x,
#    '(MeV) ',1x,'(mb/sr)      ',1x,'(mb/sr)      ',1x,
#    '(1)          ',1x,'(1) ',7x,'(1) ',7x,'(ns)          ')
#    write (1,113)
# 113format('SCATTERED particle:')
#    write (1,11)
#    write (1,111)
#    write (1,112)

for theta in numpy.linspace(amin, amax, (amax - amin) / astep): # do 5 i=1,istep          ! Step through scattering angles in lab frame.
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

    # print theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio,beta1,beta2,timdi
    print theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio,beta1,beta2,timdi
# 77   format(f6.2,3x,f6.2,3x,f6.2,3x,f6.1,4x,f6.1,1x,
#              e10.3,4x,e10.3,4x,e10.3,4x,e10.3,1x,e10.3,1x,e10.3)


if thelim < 179.0:

    print 'SCATTERED particle, SECOND solution:'

    for theta in numpy.linspace(amin, thelim, (thelim - amin) / astep):

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
        # print theta,thetb,tcm,dca,e1,ruthc,ruthl,ratio

#     12 format('theta',4x,'theta',4x,'theta',4x,'dist of  ',1x,
#        'energy',1x,'dsigma/domega',1x,'dsigma/domega',1x,
#        'dsigma/domega')
#     121format('rec  ',4x,'scat ',4x,'rec  ',4x,'closest  ',1x,
#        'in lab',1x,'c.o.m.       ',1x,'laboratory   ',1x,
#        'lab/c.o.m.   ')
#     122format('lab  ',4x,'lab  ',4x,'c.o.m.',3x,'appr (fm)',1x,
#        '(MeV) ',1x,'(mb/sr)      ',1x,'(mb/sr)      ',1x,
#        '(1)          ')

# Recoiling target.
print 'RECOILING particle:'
for theta in numpy.linspace(amin, min(amax,89.), (min(amax,89.) - amin) / astep):
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

    # print theta,thetb,tcm,dca,e2,ruthc,ruthl,ratio

    # Add to lists for delta-TOF plot:
    recoil_theta_list.append(theta)
    recoil_timdi_list.append(timdi)


for i in range(len(recoil_theta_list)):
    print '(' + str(round(recoil_theta_list[i])) + ', ' + str(round(recoil_timdi_list[i])) + ')'
for i in range(len(theta_list)):
    print '(' + str(round(theta_list[i])) + ', ' + str(round(timdi_list[i])) + ')'

the_plot = plt.plot(theta_list, timdi_list)
the_plot = plt.plot(recoil_theta_list, recoil_timdi_list)
plt.xlabel('scatter / recoil angle (deg)')
plt.ylabel('time-of-flight difference (ns)')
pylab.show()



#  Rachel/rachel.py:                plt.text(text_coordinate_x, text_coordinate_y, yield_label, dict(color=circle_color, multialignment="right", rotation=text_angle))
#  Rachel/rachel.py:                #plt.text(text_coordinate_x, text_coordinate_y, yield_label, dict(color=circle_color, multialignment="left", rotation=text_angle, rotation_mode = "anchor"))
#      Rachel/rachel.py:                plt.annotate("", xy=(float_final_band_number,final_energy),  xycoords='data', xytext=(float_initial_band_number,initial_energy),\
#              Rachel/rachel.py:                plt.text(text_coordinate_x, text_coordinate_y, yield_label, dict(color=arrow_color, multialignment="center", rotation=text_angle))
#              Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#              Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#              Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#              Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#              Rachel/rachel.py:        plt.title(title)
#              Rachel/rachel.py:        plt.xlabel(x_label)
#              Rachel/rachel.py:        plt.ylabel(y_label)
#              Rachel/rachel.py:            plt.ylim(low_y_plot_limit,high_y_plot_limit)
#              Rachel/rachel.py:            plt.ylim(high_y_plot_limit,low_y_plot_limit)
#              Rachel/rachel.py:        plt.xlim(low_x_plot_limit,high_x_plot_limit)
#              Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#              Rachel/rachel.py:        plt.plot(wxlist,wylist,'r--')
#              Rachel/rachel.py:        plt.semilogy()
#              Rachel/rachel.py:        plt.scatter(spinlist,blist,s=25,c='r',marker=(0,3,0))
#              Rachel/rachel.py:        plt.ylim(y_min,y_max)
#              Rachel/rachel.py:        plt.xlim(x_min,x_max)
#              Rachel/rachel.py:        plt.xlabel('Initial spin')
#              Rachel/rachel.py:        plt.ylabel('B(ML;up)')
#              Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#              Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#              Rachel/rachel.py:        plt.annotate(arrowlabel, xy=(x2,y2),  xycoords='data',
#                  Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#                  Rachel/rachel.py:        plt.scatter(x_position,levelenergy,s=25,c=requestedcolor,marker=(0,3,0))
#                  Rachel/rachel.py:            plt.annotate(comment, xy=(x_position,levelenergy),  xycoords='data',
#                      Rachel/rachel.py:        plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#                      Rachel/rachel.py:        plt.xticks(numpy.arange(nbands + ASMIDGE),bandlabels)
#                      Rachel/rachel.py:        plt.ylim(-maxenergy/5., 1.1*maxenergy)   # a little extra head room and room 
#                      Rachel/rachel.py:        plt.xlim(0., nbands + 1.1)  # Add a little extra to round nbands up to next 
#                      Rachel/rachel.py:            plt.figure(LEVELSCHEMEFIGURE,figsize=LSFIGSIZE)
#                      Rachel/rachel.py:        lowlim,highlim = plt.ylim()
# 
# 

