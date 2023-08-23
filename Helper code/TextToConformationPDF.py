# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 15:52:34 2021

@author: jordy
"""

'''
Converts a text-based representation of a templated copolymerisation 
conformation to a 'blob and rod' which is then saved as a pdf.
'''

import numpy as np
import matplotlib.pyplot as plt

def xy(r,phi):
  return r*np.cos(phi), r*np.sin(phi)

def text_to_pic(lines,filename):
    fig,ax = plt.subplots(dpi = 500,figsize=(10,10))
   
    offset = 6
    
    # Circle format
    r = 0.4 # Radius copies
    rs = 0.45 # width templates
    
    cbw = 0.08 # Circle border width
    col0 = (240/256,240/256,240/256) # Color monomer 0
    col1 = col0#(76/256,132/256,167/256) # Color monomer 1
    cols = [col0,col1]
    mons = ['0','1']
    bordcol = 'k' # Border color
    
    tilt = 0.5#-0.5 # Tilt to tail monomers
    tailbondscale = 1.8 # Height between tail bonds
    
    ctheight = 1.6 # Copy height
    
    bondw = 0.2 # Bond thickness
    bondcol = (60/256,60/256,60/256) # Bond color
    
    generic = True # add a leading edge generic bond and specific
    genericcol = bondcol
    speccol = 'r'
    
    endoftemplate = True # Draw colour on end of template
    endcol = (210/256,245/256,50/256)
    
    d = 0.1 # Generic specific bond spacing
    specticks = 5 # Number of spec bond ticks
    
    linetemp = []
    maxlen = max([len(l) for l in lines])
    # Find the template layer labels
    for row,line in enumerate(lines):
        linetemp.append(line+(3+maxlen-len(line))*' ')
        
        if '|' in line:
            temp_row = row 
            
    lines = linetemp
    print(linetemp)
    # Loop through the rows
    for row,line in enumerate(lines):
        
        y = temp_row - row - 1
        
        if row < temp_row -1:
            # If we are in a tail row then we need to draw the angled tails
            y = tailbondscale*y
            for col,char in enumerate(line):
                
                if char in mons:
                    x = 2*( col - offset)/3
                    
                    if (col - offset)%3 == 0:
                        x = x - tilt
                    if (col - offset)%3 == 1:
                        x = x + tilt - 2/3
                    
                    ax.add_patch(plt.Circle((x, y), r-cbw, facecolor=cols[int(char)],zorder = 100))
                    ax.add_patch(plt.Circle((x, y), r, facecolor=bordcol,zorder = 90))
                    
                    # Now add the backbones connecting tails
                    if lines[row+1][col] in mons:
                        
                        # Tilted backbone for first units in tail
                        if row == temp_row - 2:
                            # Lagging tails
                            if (col - offset)%3 == 0:
                                ax.add_patch(plt.Rectangle((x-bondw/2, y),\
                                                            width = bondw, \
                                                            height = -np.sqrt(tilt**2 + tailbondscale**2),\
                                                            angle = 180/np.pi*np.arctan((tilt-bondw/2)/tailbondscale), \
                                                            fill = True,facecolor=bondcol,zorder = -100))
                                
                                                        
                        else:
                            ax.add_patch(plt.Rectangle((x-bondw/2, y), 
                                                            width = bondw, 
                                                            height = -tailbondscale,
                                                            fill = True,facecolor=bondcol,zorder = -100))
                            
                    elif lines[row+1][col-1] in mons:
                        # Leading tail
                        if row == temp_row - 2:
                            
                            ax.add_patch(plt.Rectangle((x-bondw/2, y), 
                                                        width = bondw, 
                                                        height = -np.sqrt(tilt**2 + tailbondscale**2),
                                                        angle = -180/np.pi*np.arctan((tilt-bondw/2)/tailbondscale), 
                                                        fill = True,facecolor=bondcol,zorder = -100))
                
                
        
        if row == temp_row - 1:
            # Now we're on the copy template connected layer
            for col,char in enumerate(line):
                
                x = 2*( col - offset)/3
                
                if char in mons:
                    
                    ax.add_patch(plt.Circle((x, y), r-cbw, facecolor=cols[int(char)],zorder = 100))
                    ax.add_patch(plt.Circle((x, y), r, facecolor=bordcol,zorder = 90))
                    #    ax.add_patch(plt.Circle((x, y), r, facecolor=cols[int(char)],edgecolor = bordcol,linestyle = '-',linewidth = clw))
                                    
                if line[col:col+2] == '--':
                    ax.add_patch(plt.Rectangle((x-1, y-bondw/2), 
                                                    width = 2, 
                                                    height = bondw,
                                                    fill = True,facecolor=bondcol,zorder = -100))
                    #ax.plot([x-1+r,x+2-r],[y,y],linewidth = bondlw,color = bondcol, zorder=-100) 
            
                
        if row == temp_row :
            # Now below the copy template connected layer
            for col,char in enumerate(line):
                
                x = 2*( col - offset)/3
                if char == '|':
                    if generic == False:
                        ax.add_patch(plt.Rectangle((x-bondw/2, 0), 
                                                            width = bondw, 
                                                            height = -ctheight,
                                                            fill = True,facecolor=bondcol,zorder = -100))
                        
                        #ax.plot([x,x],[-ctheight,0],linewidth = bondlw,color = bondcol, zorder=-100)
                        #ax.add_patch(plt.Rectangle((x-bondlw/2, y-1+r-bondlw/2), bondlw,y+1-2*r+bondlw/2, facecolor=bondcol,edgecolor = 'k',linestyle = '-'))
                
                    # Draw a generic bond if appropriate
                    elif generic == True:
                        
                        # If we aren't at a leading edge, draw dotted line
                        if (lines[row-2][col+1] in mons) | (lines[row-1][col+1:col+3] == '--'):
                            
                            for i in range(specticks):
                                ax.add_patch(plt.Rectangle((x-bondw/2, -i*ctheight/specticks), 
                                                            width = bondw, 
                                                            height = -ctheight/specticks/2,
                                                            fill = True,facecolor=speccol,zorder = -100))
                                
                        else: # If we are at leading edge draw dotted and solid line displaced by a small amounnt from center
                            
                            for i in range(specticks):
                                ax.add_patch(plt.Rectangle((x-bondw/2-d, -i*ctheight/specticks), 
                                                            width = bondw-d/2, 
                                                            height = -ctheight/specticks/2,
                                                            fill = True,facecolor=speccol,zorder = -100))
                            ax.add_patch(plt.Rectangle((x+bondw/2+d, 0), 
                                                            width = -bondw+d/2, 
                                                            height = -ctheight,
                                                            fill = True,facecolor=genericcol,zorder = -100))
                            
                            # ax.plot([x-d,x-d],[-ctheight,0],':',linewidth = bondlw,color = bondcol, zorder=-100)
                            # ax.plot([x+d,x+d],[-ctheight,0],'-',linewidth = bondlw,color = bondcol, zorder=-100)
        
        if row == temp_row +1 :
            # Now Template
            for col,char in enumerate(line):
                
                x = 2*( col - offset)/3    
                
                # Bit of logic to not draw lines if at end of template
                if line[col:col+4] == '0  0':
                    
                    ax.add_patch(plt.Rectangle((x, -ctheight-bondw/2), 
                                                    width = 2, 
                                                    height = bondw,
                                                    fill = True,facecolor=bondcol,zorder = -100))
                                        
                if char in mons:
                    ax.add_patch(plt.Rectangle((x-rs+cbw, -ctheight-rs+cbw), 
                                                width = 2*rs-2*cbw,
                                                height = 2*rs-2*cbw, 
                                                fill = True,facecolor=cols[int(char)],zorder =100))
                    
                    ax.add_patch(plt.Rectangle((x-rs, -ctheight-rs), 
                                                width = 2*rs,
                                                height = 2*rs, 
                                                fill = True, facecolor=bordcol,zorder = 90))
                    # Draw end of template colour box                    
                    if (endoftemplate == True) & (line[col:col+4] != '0  0'):
                        ax.add_patch(plt.Rectangle((x-rs+cbw, -ctheight-rs+cbw), 
                                                width = 2*rs-2*cbw,
                                                height = 2*rs-2*cbw, 
                                                fill = True,facecolor=endcol,zorder =100))
                    
                        ax.add_patch(plt.Rectangle((x-rs, -ctheight-rs), 
                                                    width = 2*rs,
                                                    height = 2*rs, 
                                                    fill = True, facecolor=bordcol,zorder = 90))
                        
                        
    ax.axis('equal')
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(filename + ".pdf", format = 'pdf', dpi=500)
    
    
with open("conformations.txt") as f:
    # Import confs text files
    lines = f.readlines()
    #print(lines)
    
    labelrows = []
    filenames = []
    # Find the start and end of each conf
    for row, line in enumerate(lines):
        if 'Title ' in line:
            labelrows.append(row)
            filenames.append(line[6:-1])
    
    
    # plot each of the conformations
    for i in range(len(labelrows)-1):
        text_to_pic(lines[labelrows[i]+1:labelrows[i+1]],filenames[i])
    text_to_pic(lines[labelrows[-1]+1:],filenames[-1])

#ax.set_yticks([])
#ax.set_xticks([])