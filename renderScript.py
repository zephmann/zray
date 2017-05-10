import os, re, sys
import maya.cmds as cmds

global zRayWind
global aaCheckBox
global aaSamp
global raySamp
global minNoise
global totalDisp
global minDepth
global maxDepth
global diffDepth
global fileHeader
global fullPath
global startFrame
global endFrame
global frameRange


# create UI window ##################################################
def zRayUI():
    global zRayWind
    global aaCheckBox
    global aaSamp
    global raySamp
    global minNoise
    global totalDisp
    global minDepth
    global maxDepth
    global diffDepth
    global fileHeader
    global fullPath
    global startFrame
    global endFrame
    global frameRange
    
    # if old window exists, delete it
    if cmds.window("zray", ex=1):
        cmds.deleteUI("zray")
    
    if cmds.windowPref("zray", ex=1):
        cmds.windowPref("zray", r=1)
    
    # create new window
    zRayWind = cmds.window("zray", t="~zRay~", widthHeight=(600, 375))
    
    # create sampling frame
    cmds.frameLayout(bs="etchedIn", l="Sampling")
    cmds.columnLayout()
    
    # sliders for pixel samples
    aaCheckBox = cmds.checkBox(l="Anti-Aliased", v=1, cc='updateSamps()')
    aaSamp = cmds.intSliderGrp(l="AA Samples", f=1, min=1, max=16, v= 4, en=1, cc='updateSamps()', dc='updateSamps()')
    raySamp = cmds.intSliderGrp(l="Ray Samples", f=1, min=1, max=10000, v=16, cc='updateSamps()', dc='updateSamps()')
    minNoise = cmds.floatSliderGrp(l="Min Noise", f=1, min=0.0, max=1.0, v=0.005, pre=6)
    
    # display total number of pixel samples
    cmds.rowLayout(nc=2)
    cmds.text(l="Total Pixel Samples ")
    totalDisp = cmds.intField(ed=0, v=256)
    cmds.setParent('..')
    
    cmds.setParent('..')
    
    # create ray limit frame
    cmds.frameLayout(bs="etchedIn", l="Ray Depth")
    cmds.columnLayout()
    
    minDepth = cmds.intSliderGrp(l="Min Ray Depth", f=1, min=1, max=10, v=3)
    maxDepth = cmds.intSliderGrp(l="Max Ray Depth", f=1, min=1, max=20, v=6)
    diffDepth = cmds.intSliderGrp(l="Direct Lighting Limit", f=1, min=1, max=20, v=2)
    
    cmds.setParent('..')
    
    # create output file frame
    cmds.frameLayout(bs="etchedIn", l="File")
    cmds.columnLayout()
    
    cmds.rowLayout(nc=2)
    cmds.text(l="Filename ")
    
    # get scene name
    sceneName = os.path.splitext(os.path.basename(cmds.file(q=True, sn=True)))[0]
    if sceneName == "":
        sceneName = "untitled"
    
    name = "frame.ppm"
    fileHeader = cmds.textField(width=330, tx=name, cc="updateFullPath()")
    cmds.setParent('..')

    cmds.rowLayout(nc=2)
    cmds.text(l="Full Path: ")
    outDir = cmds.workspace(q=1, dir=1)
    fullPath = cmds.textField(ed=0, width=525, tx=outDir+"zray/"+sceneName+"/images/"+name)
    cmds.setParent('..')
    
    frameRange = cmds.checkBox(l="Frame Range", v=0, cc='updateFrameRange()')
    
    cmds.rowLayout(nc=2)
    startFrame = cmds.intField(ed=1, v=1, en=0)
    endFrame = cmds.intField(ed=1, v=240, en=0)
    cmds.setParent('..')
    
    renderBtn = cmds.button(l="Render", c='zRender()', width=200)
    
    cmds.setParent('..')
    
    cmds.showWindow(zRayWind)


def updateSamps():
    global zRayWind
    global aaCheckBox
    global aaSamp
    global raySamp
    global minNoise
    global totalDisp
    global minDepth
    global maxDepth
    global diffDepth
    global fileHeader
    global fullPath
    global startFrame
    global endFrame
    global frameRange
        
    totalSamp = 1
    activeVal = 0
    if cmds.checkBox(aaCheckBox, q=1, v=1) == 1:
        totalSamp = cmds.intSliderGrp(aaSamp, q=1, v=1)
        activeVal = 1
    totalSamp = totalSamp * totalSamp * cmds.intSliderGrp(raySamp, q=1, v=1)
    
    cmds.intSliderGrp(aaSamp, e=1, en=activeVal)
    cmds.intField(totalDisp, e=1, v=totalSamp)


def updateFrameRange():
    global zRayWind
    global aaCheckBox
    global aaSamp
    global raySamp
    global minNoise
    global totalDisp
    global minDepth
    global maxDepth
    global diffDepth
    global fileHeader
    global fullPath
    global startFrame
    global endFrame
    global frameRange
    
    activeVal = cmds.checkBox(frameRange, q=1, v=1)
    
    cmds.intField(startFrame, e=1, en=activeVal)
    cmds.intField(endFrame, e=1, en=activeVal)


def updateFullPath():
    global zRayWind
    global aaCheckBox
    global aaSamp
    global raySamp
    global minNoise
    global totalDisp
    global minDepth
    global maxDepth
    global diffDepth
    global fileHeader
    global fullPath
    global startFrame
    global endFrame
    global frameRange
    
    name = cmds.textField(fileHeader, q=1, tx=1)
    outDir = cmds.workspace(q=1, dir=1)
    sceneName = os.path.splitext(os.path.basename(cmds.file(q=True, sn=True)))[0]
    if sceneName == "":
        sceneName = "untitled"
    
    cmds.textField(fullPath, e=1, tx=outDir+"zray/"+sceneName+"/images/"+name)



def zRender():
    global zRayWind
    global aaCheckBox
    global aaSamp
    global raySamp
    global minNoise
    global totalDisp
    global minDepth
    global maxDepth
    global diffDepth
    global fileHeader
    global fullPath
    global startFrame
    global endFrame
    global frameRange
    
    fn = str(cmds.textField(fileHeader, q=1, tx=1))
    aa = int(cmds.checkBox(aaCheckBox, q=1, v=1))
    aaS = int(cmds.intSliderGrp(aaSamp, q=1, v=1))
    rayS = int(cmds.intSliderGrp(raySamp, q=1, v=1))
    noise = float(cmds.floatSliderGrp(minNoise, q=1, v=1))
    minD = int(cmds.intSliderGrp(minDepth, q=1, v=1))
    maxD = int(cmds.intSliderGrp(maxDepth, q=1, v=1))
    difD = int(cmds.intSliderGrp(diffDepth, q=1, v=1))
    
    cur = int(cmds.currentTime(q=1))
    
    finalFile = fn
    
    ######################################################################
    ######################################################################
    ## JOE, before giving to students: ###################################
    ## comment the line below ############################################
    
    #renderPath = "/home/" + os.getenv("USERNAME") + "/mount/fachome/zRay/zRay.out"
    
    ## Then uncomment out this line ########################################
    renderPath = "/home/" + os.getenv("USERNAME") + "/mount/stuhome/zRay/zRay.out"
    
    ######################################################################
    ######################################################################
    ######################################################################
    
    l = [m.start() for m in re.finditer('#', fn)]
    length = len(l)
    
    # if rendering a sequence
    if cmds.checkBox(frameRange, q=1, v=1) == 1:
        
        # warn if there are no wild cards
        if length == 0:
            print "#No Wild Cards#"
            return
        
        # else loop through frames, replace wild cards then call render function
        else:
            
            # flist and digits should be lists of the substrings, and the number of pound signs
            flist = []
            digits = []
            flength = 0
            dlength = 0
            
            # if there is only group of wild cards
            if length-1 == l[length-1] - l[0]:
                for s in fn.split('#'):
                    if not s == "":
                        flist.append(s)
                flength = 2
                digits = [l[length-1] - l[0] + 1]
                dlength = 1
            
            # else if there are multiple groups
            else:
                for s in fn.split('#'):
                    if not s == "":
                        flist.append(s)
                        flength += 1
                
                index = 0
                while index < length:
                    start = index
                    
                    while index < length-1 and l[index+1] - l[index] == 1:
                        index += 1
                    
                    digits.append(index-start+1)
                    dlength += 1
                    index += 1
            
            # loop through frame range, update file name
            start = int(cmds.intField(startFrame, q=1, v=1))
            end = int(cmds.intField(endFrame, q=1, v=1))
            
            totalFrames = abs(end-start)
            
            frame = start
            
            # if frame range is ascending
            if start < end:
                while frame <= end:
                    
                    # updated file name
                    finalFile = ""
                    
                    # if file name doesn't start with a wild card add first substring
                    count = 0
                    if not l[0] == 0:
                        finalFile = flist[0]
                        count = 1
                    
                    # loop through digits and add the frame number padded with the right
                    # number of zeroes
                    i = 0
                    while i < dlength:
                        finalFile = finalFile + str(frame).zfill(digits[i])
                        
                        # if we haven't exhausted the list of substrings, add the next one
                        if count < flength:
                            finalFile = finalFile + flist[count]
                            count += 1
                        
                        i+=1
                    
                    # call render function
                    print ("rendering: %s" % finalFile)
                    
                    # wait til render is finished, then move to next frame
                    
                    # move to next frame
                    frame += 1
            
            else:
                while frame >= end:
                     
                     # updated file name
                    finalFile = ""
                    
                    # if file name doesn't start with a wild card add first substring
                    count = 0
                    if not l[0] == 0:
                        finalFile = flist[0]
                        count = 1
                    
                    # loop through digits and add the frame number padded with the right
                    # number of zeroes
                    i = 0
                    while i < dlength:
                        finalFile = finalFile + str(frame).zfill(digits[i])
                        
                        # if we haven't exhausted the list of substrings, add the next one
                        if count < flength:
                            finalFile = finalFile + flist[count]
                            count += 1
                        
                        i+=1
                    
                    # call render function
                    print ("rendering: %s" % finalFile)
                    
                    # wait til render is finished, then move to next frame
                    
                    # move to next frame
                    frame -= 1
        
        # reset to current frame
        cmds.currentTime(cur)
    
    # else if single frame
    else:
        
        # if file contains wildcards, replace them with current frame number
        if length > 0:
            
            # flist and digits should be lists of the substrings, and the number of pound signs
            flist = []
            digits = []
            flength = 0
            dlength = 0
            
            # if there is only group of wild cards
            if length-1 == l[length-1] - l[0]:
                for s in fn.split('#'):
                    if not s == "":
                        flist.append(s)
                flength = 2
                digits = [l[length-1] - l[0] + 1]
                dlength = 1
            
            else:
                for s in fn.split('#'):
                    if not s == "":
                        flist.append(s)
                        flength += 1
                
                index = 0
                while index < length:
                    start = index
                    
                    while index < length-1 and l[index+1] - l[index] == 1:
                        index += 1
                    
                    digits.append(index-start+1)
                    dlength += 1
                    index += 1
            
            # replace wildcards with padded frame number
            finalFile = ""
            
            # if file name doesn't start with a wild card add first substring
            count = 0
            if not l[0] == 0:
                finalFile = flist[0]
                count = 1
            
            # loop through digits and add the frame number padded with the right
            # number of zeroes
            i = 0
            while i < dlength:
                finalFile = finalFile + str(cur).zfill(digits[i])
                
                # if we haven't exhausted the list of substrings, add the next one
                if count < flength:
                    finalFile = finalFile + flist[count]
                    count += 1
                
                i+=1
        
        # call render function and return
        outDir = cmds.workspace(q=1, dir=1)
        sceneName = os.path.splitext(os.path.basename(cmds.file(q=True, sn=True)))[0]
        if sceneName == "":
            sceneName = "untitled"
        fn = outDir + finalFile
        
        # generate .dor file
        print ("rendering: %s" % finalFile)
        
        dorname, imgname = zRaySceneGen(fn, minD, maxD, difD, aa, aaS, rayS, noise)
        
        # call renderer
        os.system(renderPath + " " + dorname)
        os.system("eog " + imgname + "&")
    
    
    ## End single frame ###################################
    
    print "\nFinished rendering\n"
    return
    
    outDir = cmds.workspace(q=1, dir=1)
    sceneName = os.path.splitext(os.path.basename(cmds.file(q=True, sn=True)))[0]
    if sceneName == "":
        sceneName = "untitled"
    fn = outDir + fn
    
    print 'zRaySceneGen(\"' + fn + '\", ' + str(minD) + ', ' + str(maxD) + ', ' + str(difD) + ', ' + str(aa) + ', ' + str(aaS) + ', ' + str(rayS) + ', ' + str(noise) + ')'
    print "Result: " + str(zRaySceneGen(fn, minD, maxD, difD, aa, aaS, rayS, noise))
# END UI ###################################################################################################################################


############################################################################################################################################
# Scene Description Generation ############################################################################################################
def zRaySceneGen(outname, minD, maxD, difD, antiAlias, aaSamp, raySamp, noise):
    # open scene description file
    filename = os.path.splitext(outname)[0] + ".dor"
    
    sceneName = os.path.splitext(os.path.basename(cmds.file(q=True, sn=True)))[0]
    if sceneName == "":
        sceneName = "~untitled~"
    
    fsplit = os.path.split(filename)
    zrayDir = fsplit[0] + "/zray/"
    
    
    if not os.path.exists(zrayDir):
        os.mkdir(zrayDir)
    
    dorDir = zrayDir + "dor/"
    imgDir = zrayDir + "images/"
    
    if not os.path.exists(dorDir):
        os.mkdir(dorDir)
    
    if not os.path.exists(imgDir):
        os.mkdir(imgDir)
    
    dorname = dorDir + fsplit[1]
    imgname = imgDir + os.path.split(outname)[1]    
    
    sceneFile = open(dorname, 'w')
    
    # get scene name and image dimensions
    width = cmds.getAttr('defaultResolution.width')
    height = cmds.getAttr('defaultResolution.height')
    
    #sceneFile.write("%s\n%s \n" % (sceneName, outname))
    sceneFile.write("%s\n%s \n" % (sceneName, imgname))
    sceneFile.write("%d %d\n" % (width, height))
    sceneFile.write("%d %d %d\n" % (minD, maxD, difD))
    sceneFile.write("%d %d %d\n%f\n\n" % (not antiAlias, aaSamp, raySamp, noise))
    
    curSel = cmds.ls(sl=1)
    ##############################################################################
    
    
    # get camera info
    cams = cmds.ls(ca=True)
    for c in cams:
        if cmds.getAttr(c + '.renderable'):
            camShape = c
    
    cam = str(cmds.listRelatives(camShape, p=True, f=True)[0])
    
    camPos = (cmds.getAttr(cam + '.translateX'), cmds.getAttr(cam + '.translateY'), cmds.getAttr(cam + '.translateZ'))
    camRot = (cmds.getAttr(cam + '.rotateX'), cmds.getAttr(cam + '.rotateY'), cmds.getAttr(cam + '.rotateZ'))
    angleOfView = cmds.camera(camShape, q=True, hfv=True)
    envCol = cmds.getAttr(camShape + '.backgroundColor')[0]
    
    sceneFile.write("Camera:\nPosition " + "%f "*3%camPos + "\nRotation " + "%f "*3%camRot)
    sceneFile.write("\nAoV %f\nBackground"%angleOfView + " %f"*3%envCol + "\n\n")
    ##############################################################################
    
    
    # get lights info
    lights = cmds.ls(lt=True)
    lightDict = {'directionalLight': 0, 'pointLight': 1, 'spotLight': 2}
    
    sceneFile.write("Num_Lights %d\n\n" % len(lights))
    
    for lightShape in lights:
        light = str(cmds.listRelatives(lightShape, p=True, f=True)[0])
        
        ligType = cmds.nodeType(lightShape)
        
        sceneFile.write('Name %s\nType %s\n'% (light[1:],lightDict[ligType]))
        
        ligRot = (cmds.getAttr(light + '.rotateX'), cmds.getAttr(light + '.rotateY'), cmds.getAttr(light + '.rotateZ'))
        ligColor = cmds.getAttr(lightShape + '.color')[0]
        ligInten = cmds.getAttr(lightShape + '.intensity')
        
        sceneFile.write('Rotation ' + '%f '*3%ligRot + '\nColor ' + '%f '*3%ligColor + '\nIntensity %f\n'%ligInten)
        
        ligDiff = int(cmds.getAttr(lightShape + '.emitDiffuse'))
        ligSpec = int(cmds.getAttr(lightShape + '.emitSpecular'))
        ligShad = int(cmds.getAttr(lightShape + '.useRayTraceShadows'))
        
        sceneFile.write("Diffuse %d\nSpecular %d\nShadows %d\n" % (ligDiff, ligSpec, ligShad))
        
        if ligType != 'directionalLight':
            ligRad = cmds.getAttr(lightShape + '.lightRadius')
        else:
            ligRad = cmds.getAttr(lightShape + '.lightAngle')
        shadRay = cmds.getAttr(lightShape + '.shadowRays')
        shadDep = cmds.getAttr(lightShape + '.rayDepthLimit')
           
        
        sceneFile.write("LightRadius %f\nShadowRays %d\nRayDepth %d\n" % (ligRad, shadRay, shadDep))
        
        if ligType != 'directionalLight':
            ligPos = (cmds.getAttr(light + '.translateX'), cmds.getAttr(light + '.translateY'), cmds.getAttr(light + '.translateZ'))
            ligRate = cmds.getAttr(lightShape + '.decayRate')
            
            sceneFile.write('Translation ' + '%f '*3%ligPos + '\nDecay %d\n'%ligRate)
            
            if ligType == 'spotLight':
                ligCone = cmds.getAttr(lightShape + '.coneAngle')
                ligPen = cmds.getAttr(lightShape + '.penumbraAngle')
                ligDrop = cmds.getAttr(lightShape + '.dropoff')
                
                sceneFile.write('Cone %f\nPenumbra %f\nDropOff %f\n' % (ligCone, ligPen, ligDrop))
        
        sceneFile.write('\n')
    ##############################################################################
    
    
    # get texture files
    files = cmds.ls(type = "file")
    textures = [cmds.getAttr(f + ".fileTextureName") for f in files]
    textures = list(set(textures))
    
    sceneFile.write('Num_Textures %d\n\n'%len(textures))
    
    for t in textures:
        sceneFile.write(t + '\n')
    
    ##############################################################################
    
    
    # get materials info
    mats = cmds.ls(mat=True)
    matDict = {'lambert': 0, 'phong':1}
    
    matCount = 0
    for m in mats:
        type = cmds.nodeType(m)
        if type == 'lambert' or type == 'phong':
            matCount += 1
    
    sceneFile.write("Num_Materials %d\n\n" % matCount)
    
    for m in mats:
        type = cmds.nodeType(m)
        if type == 'lambert' or type == 'phong':
            matDiff = cmds.getAttr(m + '.diffuse')
            
            sceneFile.write('Name %s\nType %d\nDiffuse %f\nColor '%(m, matDict[type], matDiff))
            
            colorConnections = cmds.ls(cmds.listConnections(m + ".color"), type = "file")
            if colorConnections:
                matColor = cmds.getAttr(colorConnections[0] + ".fileTextureName")
                texIndex = textures.index(matColor)
                sceneFile.write('Texture %d\n' % texIndex)
            else:
                matColor = cmds.getAttr(m + '.color')[0]
                sceneFile.write('%f '*3%matColor + '\n')
            
            colorConnections = cmds.ls(cmds.listConnections(m + ".transparency"), type = "file")
            if colorConnections:
                matTrans = cmds.getAttr(colorConnections[0] + ".fileTextureName")
                texIndex = textures.index(matTrans)
                sceneFile.write('Transparency Texture %s\n' % texIndex)
            else:
                matTrans = cmds.getAttr(m + '.transparency')[0]
                matTrans = (1-matTrans[0], 1-matTrans[1], 1-matTrans[2])
                sceneFile.write('Transparency ' + '%f '*3%matTrans + '\n')
            
            colorConnections = cmds.ls(cmds.listConnections(m + ".ambientColor"), type = "file")
            if colorConnections:
                matAmbi = cmds.getAttr(colorConnections[0] + ".fileTextureName")
                texIndex = textures.index(matAmbi)
                sceneFile.write('Ambient Texture %d\n' % texIndex)
            else:
                matAmbi = cmds.getAttr(m + '.ambientColor')[0]
                sceneFile.write('Ambient ' + '%f '*3%matAmbi + '\n')
            
            colorConnections = cmds.ls(cmds.listConnections(m + ".incandescence"), type = "file")
            if colorConnections:
                matIncan = cmds.getAttr(colorConnections[0] + ".fileTextureName")
                texIndex = textures.index(matIncan)
                sceneFile.write('Incandescence Texture %d\n' % texIndex)
            else:
                matIncan = cmds.getAttr(m + '.incandescence')[0]
                sceneFile.write('Incandescence ' + '%f '*3%matIncan + '\n')
            
            matRefr = int(cmds.getAttr(m + ".refractions"))
            matRefrIndex = cmds.getAttr(m + ".refractiveIndex")
            
            sceneFile.write('Refractions %d\nRefrIndex %f\n'%(matRefr, matRefrIndex))
            
            if type == 'phong':
                matCos = cmds.getAttr(m + '.cosinePower') * cmds.getAttr(m + '.chromaticAberration')
                sceneFile.write('Cosine %f\n' %matCos)
                
                colorConnections = cmds.ls(cmds.listConnections(m + ".specularColor"), type = "file")
                if colorConnections:
                    matSpec = cmds.getAttr(colorConnections[0] + ".fileTextureName")
                    texIndex = textures.index(matSpec)
                    sceneFile.write('Specular Texture %d' % texIndex)
                else:
                    matSpec = cmds.getAttr(m + '.specularColor')[0]
                    sceneFile.write('Specular ' + '%f '*3%matSpec)
                    
                matRefl = cmds.getAttr(m + '.reflectivity')
                sceneFile.write('\nReflectivity %f\n'%matRefl)
            
            bumpNode = cmds.listConnections(m + ".normalCamera")
            if bumpNode:
                bumpNode = bumpNode[0]
                bumpDepth = cmds.getAttr(bumpNode + ".bumpDepth")
                bumpValue = cmds.getAttr(bumpNode + ".bumpValue")
                if bumpValue == 1:
                    bumpText = cmds.getAttr(cmds.listConnections(bumpNode + ".bumpValue")[0] + ".fileTextureName")
                    sceneFile.write('Bump %s %f\n' % (bumpText, bumpDepth))
            
            sceneFile.write('\n')
    ##############################################################################
    
    
    # get mesh info
    shapes = cmds.ls(g=True)
    
    # freeze transformations
    cmds.select(shapes)
    cmds.makeIdentity( apply=True )
    
    sceneFile.write("Num_Shapes %d\n\n" % len(shapes))
    
    for s in shapes:
        print s
        cmds.select(s)
        tris = cmds.polyEvaluate(t=True)
        sceneFile.write('Mesh:\n%s\nTris %d\n' % (s,tris))
        
        faces = cmds.polyEvaluate(f=True)
        
        sceneFile.write('Indices ')
        for f in range(faces):
            cmds.select(s + '.f[' + str(f) + ']')
            vertStr = str(cmds.polyInfo(fv=True)).split()
                        
            vertStr = (vertStr[2], vertStr[3], vertStr[4])
            
            sceneFile.write('%s '*3%vertStr)
        
        cmds.select(s)
        verts = cmds.polyEvaluate(v=True)
        
        sceneFile.write('\nNum_Verts %d\nVertices '%verts)
        for v in xrange(verts):
            vert = tuple(cmds.pointPosition(s + '.vtx[%d]' % v, w=1))
            sceneFile.write('%f '*3%vert)
        
        sceneFile.write('\nNormals ')
                
        cmds.select(cmds.polyListComponentConversion(s, tf=1), r=1)
        faces = cmds.ls(sl=1,fl=1)
        
        print "normals"
        
        for f in faces:
            cmds.select(cmds.polyListComponentConversion(f, tvf=1), r=1)
            vFaces = cmds.ls(sl=1, fl=1)
            for vf in vFaces:
                norms = tuple(cmds.polyNormalPerVertex(vf, q=1, xyz=1))
                sceneFile.write('%f '*3%norms)
        
        
        print "uvs"
        
        if cmds.polyListComponentConversion(faces, toUV=1) != []:
            sceneFile.write('\nUVs ')
            
            for f in faces:
                cmds.select(cmds.polyListComponentConversion(f, tvf=1), r=1)
                vFaces = cmds.ls(sl=1, fl=1)
                for vf in vFaces:
                    cmds.select(cmds.polyListComponentConversion(vf, tuv=1), r=1)
                    uv = tuple(cmds.polyEditUV(q=1))
                    sceneFile.write('%f '*2%uv)
        else:
            sceneFile.close()
            print '%s has no UV data!'%s
            return 1
        
        double = cmds.getAttr(s + ".doubleSided")
        sceneFile.write('\nDouble %d'%double)
        
        primary = cmds.getAttr(s + ".primaryVisibility")
        sceneFile.write('\nVisible %d'%primary)
        
        castsShadows = cmds.getAttr(s + ".castsShadows")
        sceneFile.write('\nShadows %d'%castsShadows)
        
        sgNodes = cmds.listConnections(s, type='shadingEngine')
        objMat = cmds.listConnections(sgNodes [0] + '.surfaceShader')[0]
        
        type = cmds.nodeType(objMat)
        if type != 'lambert' and type != 'phong':
            objMat = 'lambert1'
        
        sceneFile.write('\nMaterial %s\n'%objMat)
        
        sceneFile.write('\n\n')
    
    sceneFile.close()
    return dorname, imgname
### END SCENE DESCRIPTION ##########################################



#zRayUI()
