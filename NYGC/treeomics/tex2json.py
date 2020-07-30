#!/usr/bin/env python 

import sys
import re

f = open(sys.argv[1], "r")
texinput = f.read()

texinput = texinput.replace("\\","")
texinput = texinput[slice(texinput.index("Tree"), texinput.index("end{tikzpicture}"))]
texinput = texinput.split('\n')
output = []
i = 0
while(i < len(texinput)):
    if("% Acquired" not in texinput[i] and "% Present" not in texinput[i]):
        if(i==0):
            output.append(texinput[i].strip().replace("Tree [.",'"name":"')+'","children":[')
        else:
            if ('footnotesize' in texinput[i] or 'NavyBlue' in texinput[i]):
                tmpstr = texinput[i].strip()
                tmpstr = re.sub(r'node\[below, black!60\]\{footnotesize\{\d*%\}\};',"",tmpstr)
                output.append(tmpstr
                            .replace("edge node[above, NavyBlue]{",'{"edge":')
                            .replace("footnotesize{",'"')
                            .replace("}}",'",')
                            .replace("};",'"",')
                            .replace(";",'')
                )
            elif ('.small' in texinput[i]):
                output.append(texinput[i].strip()
                            .replace("[.small{",'"name":"')
                            .replace("}",'","children":[')
                )
            elif ('black,draw,text' in texinput[i]):
                output.append(texinput[i].strip()
                            .replace("node[black,draw,text width=1.09cm,inner sep=2pt,align=center]{",'"name":"')
                            .replace("};", '"},')

                )
            elif ('% VAF of acquired mutations' in texinput[i]):
                output.append('"mean_vaf": "' + texinput[i][slice(texinput[i].index('mean')+ 5,texinput[i].index(';'))] + '",')
                output.append('"median_vaf": "' + texinput[i][slice(texinput[i].index('median')+ 7, len(texinput[i]))] + '",')
            else:
                output.append(texinput[i].strip()
                            .replace("]","]},",1)
                           )
    i = i + 1
    
output = "{" + "".join(output).replace("},]", "}]")
n = output.rindex(",")
output = output[slice(0,n)] + output[slice(n)].replace(",","",1)
print(output)