import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['Rad3', 'Kin28','Ccl1','Ssl1','Ssl2','Tfb1','Tfb2','Tfb3','Tfb4','Tfb5']

# Volume of each protein/domain
vol = {'Rad3':105153,\
        'Kin28':47234,\
        'Ccl1':45546,\
        'Ssl1':53743,\
        'Ssl2':106272,\
        'Tfb1':75055,\
        'Tfb2':63140,\
        'Tfb3':43342,\
        'Tfb4':41508,\
        'Tfb5':9588
      }

# Color of each protein/domain
col =  {'Rad3':'darkgreen',\
        'Kin28':'red',\
        'Ccl1':'orange red',\
        'Ssl1':'purple',\
        'Ssl2':'navy blue',\
        'Tfb1':'yellow',\
        'Tfb2':'cyan',\
        'Tfb3':'orange',\
        'Tfb4':'salmon',\
        'Tfb5':'dim gray'
        }

totvol = 0.0
for p in prots:
    totvol += vol[p]

runCommand('set bgcolor white')
i=0

#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' encloseVolume '+str(vol[p]/0.7)+' step 1') #TODO play around with threshold for individual proteins
    #runCommand('volume #'+str(i)+' transparency 0.5')
    runCommand('color '+col[p]+' #'+str(i))
    #runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str( 0.7 - i / 20.0)) 
    i += 1
    
