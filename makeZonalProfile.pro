pro makeZonalProfile, tpfilelist, ufile=ufile, display=display, $
                      primefile=primefile, bigfile=bigfile, histfile = histfile, $
                      remfile=remfile, noOvals=noOvals, orientFile=orientFile, $
                      prodFile=prodFile, centerFile=centerFile, elimfile=elimfile, $
                      longlimits=longlimits, savefile=savefile, nomask=nomask, tworots=tworots,$
                      reverse=reverse, binRight=binRight, deltay=deltay, allv=allv, vcorrect=vcorrect,$
                      ally=ally
                      

;*********************************************      
;Setup
;*********************************************
if(n_elements(vcorrect) eq 0) then vcorrect=0      
if(n_elements(display) eq 0) then display=1
if(n_elements(noOvals) eq 0) then begin   
    noOvals=0
    print, 'Ovals are being left in.'
    endif
if(n_elements(tworots) eq 0) then tworots=0
if(n_elements(nomask) eq 0) then nomask=0
if(n_elements(reverse) eq 0) then reverse=0
if(n_elements(binRight) eq 0) then binRight=0
RJ=double(7.1492e7);
Rp=double(6.6854e7);
scale = 0.1 ; degrees/pix
;linebins = indgen(90)*20+2
linebins = indgen(180)*10+2  ;runa, rung, runh, rune, runo, run s-o
boxsize=31
;linebins = indgen(180)*10+8 ;run tr-o, run ts-o
;linebins = indgen(180)*10+7 ;runq
;linebins = indgen(180)*10+9  ;runc
;linebins = indgen(90)*20+12    ;runb
;linebins = indgen(60)*30+4    ;rund
;linebins = indgen(180)*10+7    ;runf
;linebins = indgen(300)*6
totkeep=0

varZeros=0
nsigbig=4.
nlb=n_elements(linebins)
du = fltarr(nlb)
dv = fltarr(nlb)
avgu = fltarr(nlb)
su=fltarr(nlb)
sv=fltarr(nlb)
avgv = fltarr(nlb)
angmeans = fltarr(nlb)
ngsigs = fltarr(nlb)
angmeans2 = fltarr(nlb)
angsigs2 = fltarr(nlb)
dudy = fltarr(nlb)
center = fltarr(nlb)
center2 = fltarr(nlb)
skew=fltarr(nlb)
uprimevprimebar = fltarr(nlb)
sigupvpbar = fltarr(nlb)
sigu=fltarr(nlb)
sigv=fltarr(nlb)
nlist=fltarr(nlb)
smallupvpbar = fltarr(nlb)
percentElim = fltarr(nlb)
r = fltarr(nlb)
imgfile = '/home/csalyk/Jupiter/images.lst'
ntp=n_elements(tpfilelist)
nmask=makemask('n')
smask=makemask('s')
longbins = fltarr(360)
;
;*********************************************      
; Read zonal profile, list of images
;*********************************************      
readcol, '~/jplcas/vicar/calibration/zvp_jupiter.dat', temp, zvplats, zvpu, rms, uerror, $
format='d,d,d,d,d', /silent
readcol, '~/Jupiter/misc/diffs.txt', bla, longdiff, bla, timediff, /silent ;W long
if(tworots eq 1) then readcol, '~/Jupiter/misc/diffs2.txt', longdiff, timediff, /silent
;
;******for shear test only
;longdiff=dblarr(58)
;timediff=dblarr(58)+34074
;readcol, '~/Jupiter/rune/ufilee', zvplats, zvpu, zvpv, zvpdudy
;
;******
readcol, imgfile, num, name, temp, temp, date, time, type, longcenter, $
  format='i,a,a,a,a,a,a,f',  /silent
readcol, '~/Jupiter/allMaps.txt', allmaps, format='a,f', /silent
;;
;*********************************************      
;Make lists containing tiepoint pairs for ALL images.  Also, expand
;longitude and time offset lists
;*********************************************      
for i=0, ntp-1 do begin 
    print, 'reading : ', tpfilelist[i]
    readcol, tpfilelist[i], linel, sampl, liner, sampr, corr,$
      format='f,f,f,f,f', /silent
if(reverse) then readcol, tpfilelist[i], liner, sampr, linel, sampl, corr,$
      format='f,f,f,f,f', /silent

    print, i, mean(linel-liner), mean(sampl-sampr)
    filenum = getfilenum(tpfilelist[i])
    print, 'filenum = ', filenum
    if(type[filenum-1] eq 't') then mask=nmask else mask=smask
    print, 'type = ', type[filenum-1]
    
    if(noOvals eq 1) then begin
        print, 'Using oval mask'
        ovalmask=makeIndMask(allmaps[filenum-1],boxsize)
        mask=mask*ovalmask
        ;vicread, allmaps[filenum-1], data, /ash
        ;plotimage, bytscl(data*mask)
    endif

    longall = sampl*scale + longcenter[filenum-1] - 900*scale
    w=where(longall ge 360., count)
    if(count gt 0) then longall[w]-=360.
    w=where(longall le 0., count)
    if(count gt 0) then longall[w]+=360.

    if (nomask eq 1) then mask[*]=1.

;get rid of points outside of desired range
    if(n_elements(longlimits) eq 0) then begin
        keep = where(mask[round(sampl), round(linel)] * mask[round(sampr),round(liner)] $
                     eq 1, count)
    endif else begin
        keep = where(mask[round(sampl), round(linel)] * mask[round(sampr),round(liner)] $
                     eq 1 and (longall ge longlimits[0]) and (longall le longlimits[1]), count)
    endelse

print, filenum, n_elements(linel), count

    if(count gt 0) then begin
        totkeep = totkeep+count
        linel = linel[keep]
        sampl = sampl[keep]
        liner = liner[keep]
        sampr = sampr[keep]
        templongdiff = fltarr(n_elements(linel))+longdiff[filenum-1]
        temptimediff = fltarr(n_elements(linel))+timediff[filenum-1]
        tempname = strarr(n_elements(linel)) + name[filenum-1]
        templongcenterall = fltarr(n_elements(linel))+longcenter[filenum-1]

        append, linelall, linel
        append, samplall, sampl
        append, linerall, liner
        append, samprall, sampr
        append, longdiffall, templongdiff
        append, timediffall, temptimediff
        append, nameall, tempname
;        append, longcenterall, templongcenterall
    endif
endfor
;
xlengths = samprall - samplall
ylengths = linerall - linelall

;
;*********************************************      
;Convert vectors to m/s and average over longitude to get zonal flow
;*********************************************      

latbins = -linebins*scale+90.1
latall = -linelall*scale+90.1
print, mean(longdiffall)
deltalong=(xlengths/10.-longdiffall)*!pi/180.
deltalat=-ylengths/10.*!pi/180. ;****negative sign appears because N is down****

rfromcenter=RJ*Rp/sqrt(RJ^2*sin(latall*!pi/180.)^2 + Rp^2*cos(latall*!pi/180.)^2)
raxialall=rfromcenter*cos(latall*!pi/180.)

u=deltalong*raxialall/timediffall ; here deltalong is already in radians
v=deltalat*rfromcenter/timediffall ; here deltalat is already in radians

if(vcorrect eq 1) then begin  ;This is a correction, assuming all vbar is due to navigation error
    allv=dblarr(n_elements(tpfilelist))
    ally=dblarr(n_elements(tpfilelist))
    index=0
    for i=0, n_elements(tpfilelist)-1 do begin
        w=where(nameall eq nameall[index])
        print, 'mean(v), before = ', mean(v[w])
        allv[i]=mean(v[w])
        ally[i]=mean(rfromcenter[w]*deltalat[w])
        v[w]=v[w]-mean(v[w])
        index=max(w)+1
        print, 'index = ', index
        print, 'mean(v), after = ', mean(v[w])
    endfor
endif
;
;*********************************************      
;Loop through all latitude bins
;*********************************************      
for i=0, nlb-1 do begin
absdiff=abs(linerall-linebins[i])
deltay=mean(linerall-linelall)
;stop
w1=where(absdiff lt 5, count)
w2=where(linelall eq linebins[i], count)
wtemp=where(w1 ne w2, count)
if(count gt 0) then stop

if (binRight eq 1) then w=where(absdiff lt 5, count) $
else  w = where(linelall eq linebins[i], count)
   if(reverse) then w=where(linerall eq linebins[i], count)
   if (count gt 0) then begin       
        localu = u[w]
        localv = v[w]
        localname = nameall[w]
        locallinel = linelall[w]
        localsampl = samplall[w]
        localliner = linerall[w]
        localsampr = samprall[w]
        localdiffx=localsampr-localsampl
;****************************************************      
;Calculate average winds; Calculate u'v' values
;**************************************************** 
        nlist[i]=n_elements(localu)
        avgu[i] = mean(localu)
        avgv[i] = mean(localv)
        sigu[i] = stddev(localu)/sqrt(n_elements(localu)) 
        sigv[i] = stddev(localv)/sqrt(n_elements(localv)) 
        localuprime = localu - avgu[i]
        localvprime = localv - avgv[i]
        du[i] = sqrt(mean(localuprime^2.))
        dv[i] = sqrt(mean(localvprime^2.))
        uprimevprimebar[i] = mean(localuprime*localvprime)
        if(n_elements(localuprime) ge 2) then begin
            sigupvpbar[i]=stddev(localuprime*localvprime)/sqrt(n_elements(localuprime)) 
            print, '***'
            print, mean(localuprime^2)*mean(localvprime^2)
            print, mean(localuprime^2*localvprime^2)-uprimevprimebar[i]^2.+covar(localu, localv)^2.-$
              covar(localuprime^2, localvprime^2.)
            print, covar(localu, localv)^2.
            print, covar(localuprime^2., localvprime^2.)
            print, stddev(localuprime*localvprime)^2.
        endif else begin
            print, 'Variance set to zero at latitude:', latbins[i]
            sigupvpbar[i]=0
            varZeros++
        endelse
        r[i] = uprimevprimebar[i]/(du[i]*dv[i])
;****************************************************
;Orientation of u', v'
;****************************************************
        if(KEYWORD_SET(orientFile)) then begin
            localx = localuprime
            localy = localvprime
            theta = atan(localy/localx)
            w1 = where(localy lt 0 and localx lt 0)
            w2 = where(localy gt 0 and localx lt 0)
            w3 = where(localy lt 0 and localx gt 0)
            theta[w1]=!pi+theta[w1]
            theta[w2]=!pi+theta[w2]
            theta[w3]=2*!pi+theta[w3]
; Find angles of vectors around 180 degrees
            h=histogram(theta*180./!pi, locations=locations, binsize=1)
            w4=where(abs(locations-180) le 100)
            g = gaussfit(locations[w4]+.5, h[w4], a, nterms=4, estimates=[15, 180, 15, 2])
            angmeans[i]=a[1]
            angsigs[i]=a[2]/sqrt(n_elements(w4)-1)
; Find angles of vectors around 0 degrees
            w5=where(theta*180./!pi gt 180)
            theta[w5]=theta[w5]-2*!pi
            h=histogram(theta*180./!pi, locations=locations, binsize=1)
            w6=where(abs(locations) le 100)
            g = gaussfit(locations[w6]+.5, h[w6], a, nterms=4, estimates=[15, 0, 15, 2])
            angmeans2[i]=a[1]
            angsigs2[i]=a[2]/sqrt(n_elements(w6)-1)
        endif
;****************************************************
;Save all upvp values to file
;****************************************************
        if(KEYWORD_SET(savefile)) then begin
           localupvp = localuprime*localvprime
            openw, 1, savefile, /append    
            for q=0,n_elements(localupvp)-1 do begin
                printf, 1, localname[q], locallinel[q], localsampl[q], localliner[q], localsampr[q],$
                  localupvp[q], format='(a15," ", d8.3," ",d8.3," ",d8.3," ",d8.3," ",d8.3)'
            endfor
            close, 1
        endif
;****************************************************
;Find locations where uprimevprimebar is high
;****************************************************
        if(KEYWORD_SET(bigfile)) then begin
            localupvp = localuprime*localvprime
            mu = mean(localupvp)
            sd = stddev(localupvp)
            wbig = where(abs(localupvp-mu) gt nsigbig*sd, count)        
            if(count gt 0) then begin
                plothist, localupvp, binsize=.1
                oplot, [mu+nsigbig*sd, mu+nsigbig*sd], [0,100], color=1000
                oplot, [mu-nsigbig*sd, mu-nsigbig*sd], [0,100], color=1000
                cursor, x, y, wait=3                
                print, 'n_elements(wbig)', n_elements(wbig)
                upvpbig = localupvp[wbig]
                linelbig = locallinel[wbig]
                samplbig = localsampl[wbig]
                linerbig = localliner[wbig]
                samprbig = localsampr[wbig]
                namebig  = localname[wbig]
;     *********************
;     Output
;     *********************
                openw, 1, bigfile, /append    
                for k=0,n_elements(namebig)-1 do begin
                    printf, 1, namebig[k], linelbig[k], samplbig[k], linerbig[k], $
                      samprbig[k], upvpbig[k], format='(a15," ", d8.3," ",d8.3," ",d8.3," ",d8.3," ",d8.3)'
                endfor
                close, 1
            endif
        endif
;****************************************************
;Eliminate high upvp values and recalculate
;****************************************************
        if(KEYWORD_SET(elimfile)) then begin
            nsig=3.
            localupvp = localuprime*localvprime
            plothist, localupvp, xhist, yhist,binsize=.1 , /noplot
            g = gaussfit(xhist, yhist, a, nterms=4)
            sd = stddev(localupvp)
            mu = mean(localupvp)
;              oplot, [mu+nsig*sd, mu+nsig*sd], [0,100], color=1000
;              oplot, [mu-nsig*sd, mu-nsig*sd], [0,100], color=1000
;              oplot, [mean(localupvp), mean(localupvp)], [0,1000],
;              color=1000, lines=1
;              oplot, [nsig*sd, nsig*sd], [0,100], color=1000
;              oplot, [-nsig*sd, -nsig*sd], [0,100], color=1000
;              oplot, [mean(localupvp), mean(localupvp)], [0,1000], color=1000, lines=1
;              cursor, x, y, wait=3     

;            wsmall = where(abs(localupvp-mu) lt nsig*sd, count)
            wsmall = where(abs(localupvp) lt nsig*sd, count)

            percentElim[i]=1.- double(count)/double(n_elements(localupvp))
            if (count gt 0) then begin
                smallupvpbar[i] = mean(localupvp[wsmall]) 
            endif else begin
                print, 'small upvpbar set to 0'
                smallupvpbar[i]=0
            endelse
        endif
;*****************************************************
;Output orientations
;*****************************************************
        if(KEYWORD_SET(orientFile)) then begin
            openw, 1, orientFile   
            for k=0,nlb-1 do begin
                printf, 1, latbins[k], angmeans[k], angsigs[k], angmeans2[k], angsigs2[k]
            endfor
            close, 1
        endif
;*****************************************************
;Output eliminated vectors
;*****************************************************
        if(KEYWORD_SET(remfile)) then begin
            openw, 1, remfile, /append   
            for k=0,n_elements(namerem)-1 do begin
                printf, 1, namerem[k], linelrem[k], samplrem[k], linerrem[k], $
                  samprrem[k]
            endfor
            close, 1
        endif
;******************************************************
;Output histogram of u'v' values; Det u'v'bar by center of histogram
;******************************************************
;        h=histogram(localuprime*localvprime, locations=locations, binsize=.1)
;        plot, locations, h, xtitle='upvp', ytitle='Number of vectors', ps=1
;        oplot, [uprimevprimebar[i],uprimevprimebar[i]], [0,100], color=1000
;        g = gaussfit(locations+.05, h, a, nterms=4)
;        oplot, locations+.05, g, color=1000
;        center[i] = a[1]
;        sd=a[2]
;        cursor, x, y, wait=3
;        print, 'sd from fit', sd
;        print, 'const = ', a[3]
;        oplot, [center[i],center[i]], [0,100], color=1000, lines=2
;        cursor, x, y, wait=3
;        print, 'i = ', i
;        print, 'lat = ', latbins[i]
;        cursor, x, y, wait=3
;        print, ' center upvp', center[i]
;        if (mean(localuprime*localvprime) ge 20) then begin
;            plot, locations, h
;            oplot, locations+.05, g, color=1000
;            oplot, [center[i], center[i]], [0,100], color=1000
;        endif
;********************************************
;Output data from some particular latitudes (*see subroutine for details*)
;********************************************
        xvals=localuprime*localvprime
        MZP_out_for_hists, i, xvals
;********************************************
endif                 ; end of 'if there are vectors at this latitude'
endfor                          ; end of loop through latitudes
;********************************************
;Calculate du/dy        
;********************************************
for i=1, nlb-2 do begin
    dlat = latbins[i+1] - latbins[i-1]
    rfromcenter=RJ*Rp/sqrt(RJ^2*sin(latbins[i]*!pi/180.)^2 + Rp^2*cos(latbins[i]*!pi/180.)^2)
    deltay = rfromcenter*dlat*!Pi/180.
    print, rfromcenter*1.d-7
    dudy[i] = (avgu[i+1] - avgu[i-1])/deltay
endfor
;********************************************
print, 'Total number of vectors :', totkeep
print, 'Variance set to zero at ', varZeros, ' locations'
;********************************************
;Output results
;********************************************
if(display eq 1) then begin
    plot, zvplats, zvpu, ps=0, color=1000
    oplot, latbins, avgu, ps=2
    oplot, latbins, avgu, ps=0
endif

if(KEYWORD_SET(ufile)) then begin
    openw, 1, ufile
    for i=0,n_elements(latbins)-1 do begin
        printf, 1, latbins[i], ' ', avgu[i], avgv[i], dudy[i], sigu[i], sigv[i]
    endfor
    close, 1
endif

if(KEYWORD_SET(primefile)) then begin
    openw, 1, primefile    
    for i=0,n_elements(latbins)-1 do begin
        printf, 1, latbins[i], du[i], dv[i], uprimevprimebar[i], sigupvpbar[i], nlist[i], r[i], $
          format='(d8.1,d8.3,d8.3,d8.3,d8.3,i7,d8.3)'

    endfor
    close, 1
endif

if(KEYWORD_SET(elimfile)) then begin
     openw, 1, elimfile
         for i=0,n_elements(latbins)-1 do begin
        printf, 1, latbins[i], ' ',  smallupvpbar[i], percentElim[i]
    endfor
    close, 1
endif

if(KEYWORD_SET(centerFile)) then begin
    openw, 1, centerfile    
    for i=0,n_elements(latbins)-1 do begin
        printf, 1, latbins[i], ' ', du[i], dv[i], center[i], sigupvpbar[i]
    endfor
    close, 1
endif
;********************************************
;
end
