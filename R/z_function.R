# functions and test scripts JS20160713
#####################################
# Integration pixel value of region of interests
#####################################
integ.profile <- function(x, axis="H", h=c(20, 50), v=c(30, 120), disp=FALSE) {
if (disp == TRUE ){
	if (axis == "H"  ){
	ax <- 1
	plot(apply(x[h[1]:h[2], v[1]:v[2]], ax, sum),ty="l",xlab="horizontal / radial distance (pixcel)", ylab="intensity (a.u.)")
	out <- apply(x[h[1]:h[2], v[1]:v[2]], ax, sum)
	return(out)
}	else if (axis == "V" ) {
	ax <- 2
	plot(apply(x[h[1]:h[2], v[1]:v[2]], ax, sum),ty="l",xlab="vertical / azimuthal angle (degree)", ylab="intensity (a.u.)")
	out <- apply(x[h[1]:h[2], v[1]:v[2]], ax, sum)
	return(out)
}
} else {
	if (axis == "H" ){
	ax <- 1
    out <- apply(x[h[1]:h[2], v[1]:v[2]], ax, sum)
}	else if (axis == "V" ) {
	ax <- 2
	out <- apply(x[h[1]:h[2], v[1]:v[2]], ax, sum)
}	
return(out)
}
}

###############################
# noise filter (median / mean / gaussian)
###############################
noise.filter <- function(x, n=3, method="median") {
if (length(dim(x))>2){warning("data must be grayscale image")
} else{
	if(max(x)<=1){
	x <- x*(2^attr(x,"bits.per.sample")-1)
   }
	px <- nrow(x)
	py <- ncol(x)
	img <- x
	f <- n %/% 2
	x_edge_plus <- matrix(NA, px+n-1, py+n-1)
	x_edge_plus[(f+1):(f+px),(f+1):(f+py)] <- x
	rasX <- matrix(NA,length(img),(n*n))
	if (method == "median" | method == "mean") {
	for(i in 1:n){
		for(j in 1:n){
			rasX[,(n*(j-1)+i)] <- array(x_edge_plus[(i:(px+i-1)),(j:(py+j-1))])
		}
		}
		if(method == "median"){
			o.img <- apply(rasX,1,median,na.rm=T)
			dim(o.img)<- c(nrow(img),ncol(img))
		} else if(method == "mean"){
			o.img <- apply(rasX,1,mean,na.rm=T)
			dim(o.img)<- c(nrow(img),ncol(img))
		}
	} else if (method == "gaussian"){
		if(n==3){
			gcf<- c(1,2,1,2,4,2,1,2,1)/16
		} else if(n==5){
			gcf<- c(1,4,6,4,1,4,16,24,16,4,6,24,36,24,6,4,16,24,16,4,1,4,6,4,1)/256
		} else {warning("only 3 or 5 for n")
		}
		rasX <- matrix(NA,length(img),(n*n))
		for(i in 1:n){
		for(j in 1:n){
			rasX[,(n*(j-1)+i)] <- gcf[n*(j-1)+i]*array(x_edge_plus[(i:(px+i-1)),(j:(py+j-1))])
		}
		}
	o.img <- apply(rasX,1,sum,na.rm=T)
	dim(o.img)<- c(nrow(img),ncol(img))
		if(n==3){
			o.img[c(1,nrow(o.img)),c(1,ncol(o.img))] <- o.img[c(1,nrow(o.img)),c(1,ncol(o.img))]*16/9
			o.img[c(1,nrow(o.img)),2:(ncol(o.img)-1)] <- o.img[c(1,nrow(o.img)),2:(ncol(o.img)-1)]*16/12
			o.img[2:(nrow(o.img)-1),c(1,ncol(o.img))] <- o.img[2:(nrow(o.img)-1),c(1,ncol(o.img))]*16/12
		}else if(n==5){
			o.img[c(1,nrow(o.img)),c(1,ncol(o.img))] <- o.img[c(1,nrow(o.img)),c(1,ncol(o.img))]*256/121
			o.img[c(1,nrow(o.img)),c(2,ncol(o.img)-1)] <- o.img[c(1,nrow(o.img)),c(2,ncol(o.img)-1)]*256/165
			o.img[c(2,nrow(o.img)-1),c(1,ncol(o.img))] <- o.img[c(2,nrow(o.img)-1),c(1,ncol(o.img))]*256/165
			o.img[c(1,nrow(o.img)),3:(ncol(o.img)-2)] <- o.img[c(1,nrow(o.img)),3:(ncol(o.img)-2)]*256/176
			o.img[3:(nrow(o.img)-2),c(1,col(o.img))] <- o.img[3:(nrow(o.img)-2),c(1,col(o.img))]*256/176
			o.img[c(2,(nrow(o.img)-1)),c(2,(ncol(o.img)-1))] <- o.img[c(2,(nrow(o.img)-1)),c(2,(ncol(o.img)-1))]*256/189
			o.img[c(2,nrow(o.img)-1),3:(ncol(o.img)-2)] <- o.img[c(2,nrow(o.img)-1),3:(ncol(o.img)-2)]*256/240
			o.img[3:(nrow(o.img)-2),c(2,ncol(o.img)-1)] <- o.img[3:(nrow(o.img)-2),c(2,ncol(o.img)-1)]*256/240
		}
	} else {
		print(NULL)
	}
	attr(o.img, "bits.per.sample")<-attr(x,"bits.per.sample")
	attr(o.img, "samples.per.pixel")<-attr(x,"samples.per.pixel")
	out <- round(o.img)
}}
############################
# binary number to decimal #
############################
bin2dec <- function(x) {
	sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}
############################
# decimal number to binary #
############################
dec2bin <- function(x, digit=8){
  if(x <= 0 && digit <= 0){
    return(NULL)
  }else{
    return(append(Recall(x%/%2,digit-1), x%%2))
  }
}
##################################################
# counts up/down in a binary sequence (for LBP)  #
##################################################
lbnum <- function(seq) {   # routine for LBP calculation
seq.trans<- seq[2:length(seq)]
p <- seq [1:length(seq)]
q <- seq.trans[1:length(seq.trans)]
q <- c(q,seq[1])
r <- p - q
sum( r != 0)
}
############################################# 
# LBP(gray scale or binary image, r=1 or 2) #
############################################# 
lbp <- function(x, r=1) {            # main for  LBP calculation
	img <- x
	if (r!=1 && r!=2){
		return(NULL)
	} else {
		if (r<=1){
			trans <- c(-1, 0, -1,1, 0,1,1,1,1,0,1,-1,0,-1,-1,-1) ; r <- 1
		} else {
			trans <- c(-2, 1, -1,2, 1,2,2,1,2,-1,1,-2,-1,-2,-2,-1) ; r <- 2
		}
		trans <-matrix(trans,2,8)
	}
	x <- nrow(img)
	y <- ncol(img)
	o.img <- list()
	t.img <- list()
	o.img.1 <- img[(r+1):(x-r), (r+1):(y-r)]
	for ( i in 1:ncol(trans)) {
				dx <- trans[1,i]
				dy <- trans[2,i]
				o.img[[i]] <- img[(r+1+dx):(x-r+dx), (r+1+dy):(y-r+dy)] - o.img.1	
	}
ulst.d <- unlist(o.img)
ulst.d <- ifelse(ulst.d>=0, 1, 0)
bin.data <- matrix(ulst.d,length(ulst.d)/8,8)
o.img.f <- o.img.f2 <-matrix( apply(bin.data,1,bin2dec),x-2*r,y-2*r)

bn <- matrix(dec2bin(0:255),ncol=8)
tr <- apply(bn,1 ,lbnum)
tr1 <- which(tr<=2)-1
tr2 <- which(tr>2)-1

for(i in 1:58){
	o.img.f2[o.img.f2==tr1[i]] <- i-1
}
o.img.fd <- data.frame(no=1:length(o.img.f),lbp=array(o.img.f))
sub2 <- subset(o.img.fd,o.img.fd$lbp %in% tr2)
o.img.f2[sub2$no] <- 58
dim(o.img.f2) <- c(nrow(o.img.f),ncol(o.img.f))

return(list(lbp.u2 = o.img.f2, lbp.ori =o.img.f))
}
##################################################### 
#  Higher order local Autocorrelation (hlac)       #
#####################################################
hlac <- function(x, r=1, disp=FALSE) {
img <- x
nth <- c("0th","1st","2nd","3rd","4th", "5th", "6th", "7th", "8th")

i.img <- img/max(img)
dec_num <- 1:256
hlac_n <- list()
omit.bin  <- c( "00001000", "00000100", "00000010", "00000001", "11000000", "01100000", "00110000", "00011000", "00001100", "00000110", "00000011", "10000001","01011000", "00010110", "00110100", "11010000", "10000011", "00000111", "00001110", "00011100", "00111000", "01110000", "11100000", "11110000", "01111000", "00111100", "00011110", "00001111", "11000011", "10000111", "11100001", "00011111", "01111100")
o.b.n <- array()
for(i in 1:length(omit.bin)) {
	 o.b.n[i] <- bin2dec(omit.bin[i]) 
}
dec_num  <- dec_num[-o.b.n]
bin_num <- dec2bin(dec_num ,8)
hlac_all <-matrix(bin_num,length(bin_num)/8,8)
dimension <-apply(hlac_all ,1,sum)
hlac_n[[1]] <- hlac_all[which(dimension==0),]
for ( k in 1:8) {
	hlac_n[[k+1]] <- hlac_all[which(dimension==k),]
}
#
trans <- c(-r,-r,-r, 0, -r,r, 0,-r,0,r,r,-r,r,0,r,r)  #(3x3 matrix centered at (0,0))
trans <- trans + r # recentered at left-top 
trans <- matrix(trans,2,8)
bin.dat <- array()
img.bank <- list()
ini.img <- i.img[(r+1):(nrow(i.img)-r),(r+1):(ncol(i.img)-r)]  # initial image
for (i in 1:ncol(trans)){
	img.bank[[i]] <- i.img[ (trans[1,i]+1):(trans[1,i]+nrow(ini.img)),(trans[2,i]+1):(trans[2,i]+ncol(ini.img))]
}
hlac_img <- matrix(0, 223, nrow(ini.img)*ncol(ini.img))
hlac_hst <- list()
cnt <- 1
for ( dms in 1:9) {
	len <- length(hlac_n[[dms]])/8
	if ( dms ==1)   	{
		hlac <- ini.img
		hlac_hst[[dms]] <-  sum(hlac)  ##
		hlac_img[1,] <- hlac  ############ image
	} else if ( dms == 9) {
		hlac <- ini.img*img.bank[[2]]*img.bank[[3]]*img.bank[[4]]*img.bank[[1]]*img.bank[[6]]*img.bank[[7]]*img.bank[[8]]*img.bank[[5]]
		hlac_hst[[dms]]  <-  sum(hlac) ##
		hlac_img[223,] <- hlac  ############ image     		
	} else {
		hlac_out <- matrix(0, nrow(ini.img),ncol(ini.img))
		hst <- array()		
		for (m in 1:len) {
			p <- 	which(hlac_n[[dms]][m,]==1) 
			cyc <- length(p)
			cnt <- cnt + 1
			hlac <- ini.img
				for (k in 1:cyc)	{
					hlac <- hlac*img.bank[[p[k]]]
				} 
			hst[[m]] <-  sum(hlac)
			hlac_img[cnt,] <- hlac ############ image
		}
		hlac_hst[[dms]] <- hst	
	}
}
names(hlac_hst) <- nth
if (disp==TRUE) {
	return(list(features=hlac_hst,mat=hlac_img,row = nrow(ini.img), col = ncol(ini.img)))	
} else {
	return(features=hlac_hst)
}
}
##########################################################
# Gray image to binary image #
##########################################################
gray2bin <- function(x, auto=TRUE, th=200, his=FALSE, dis=FALSE) {
if (length(dim(x))>2){warning("data must be grayscale image")
} else{
if(max(x)<=1){
	x <- x*(2^attr(x,"bits.per.sample")-1)
}
if (his == TRUE){
	hist(x,main="histogram of input image",xlab="gray scale", )
	}
img.mod <- matrix(0,nrow(x),ncol(x))
if(auto){
judge<- array(0,ceiling(max(x)))
for (i in 1:ceiling(max(x))){
	wh <- x[x>=i]
	bl <- x[x<i]
	omega1 <- as.numeric(length(bl))
	m1 <- mean(bl)
	omega2 <- as.numeric(length(wh))
	m2 <- mean(wh)
	judge[i] <- omega1*omega2*abs(m1-m2)
}
judge[is.na(judge)]<-0
t <- which(judge==max(judge))
} else {
	t <- th
}
if (dis == TRUE){		
		plot(1:length(judge), judge, main="discriminant criteria", xlab="gray scale", ylab="judge",ty="l")
		abline(v=t)
	}
img.mod[x>=t] <- 1
img.mod[x<t] <- 0
attr(img.mod, "bits.per.sample")<-attr(x,"bits.per.sample")
attr(img.mod, "samples.per.pixel")<-attr(x,"samples.per.pixel")
return(img.mod)
}
}
################################ 
# RGB to Grey scale conversion #
################################
rgb2gray <- function(x, coefs=c(0.30, 0.59, 0.11)) {
	if (is.null(dim(x))) stop("image must have rgb type")
	if (length(dim(x))<3) stop("image must have rgb type")
	if (max(x) <= 1) {
		x <- x*(2^(attr(x,"bits.per.sample"))-1)
	}
	imgdata <- round((coefs[1] * x[,,1] + coefs[2] * x[,,2] + coefs[3] * x[,,3]))
	attr(imgdata, "bits.per.sample")<-attr(x,"bits.per.sample")
	attr(imgdata, "samples.per.pixel") <- 1
	return(imgdata)
}
################################### 
# Gray level cooccurence matrices  #
###################################
glcm<- function(x, t.level=4, d=1) { # GLCM 
if (length(dim(x))>2){warning("data must be grayscale image")
} else{
if(max(x)<=1){
	x <- x*(2^attr(x,"bits.per.sample")-1)
}
i.img <- x
img.info <- attributes(x)
s.level <- img.info$bits.per.sample
tgl <- 2^t.level  # gray level of target image
sgl <- 2^s.level # gray level of source image
ini.img <- matrix(round(i.img*(tgl-1)/(sgl-1)), nrow(i.img),ncol(i.img)) # linear re-scaling
#  translation in 4 directions, (0, 45, 90, 135 direction: vector length = d)
		trans <- c(d, 2*d,2*d, 0, 0, d, 0, 0)   
		trans <- matrix(trans,2,4)
		img.bank <- list()
		ref.img <- ini.img[(d+1):(nrow(ini.img)-d),(d+1):(ncol(ini.img)-d)]  # initial ROI that overlaps with translated image banks.
		o.glcm <- list()
for (i in 1:ncol(trans)){ # image bank
			img.bank[[i]] <- ini.img[ (trans[1,i]+1):(trans[1,i]+nrow(ref.img)),(trans[2,i]+1):(trans[2,i]+ncol(ref.img))]
	 }
	 	for (k in 1:ncol(trans)){
	 		ref1 <- factor(ref.img,levels=0:(tgl-1))  # i, j pair
	 		ref2 <- factor(img.bank[[k]],levels=0:(tgl-1))  # j, i pair
			t.glcm <- table(ref1,ref2) + t(table(ref1,ref2))
			m.glcm <- matrix(t.glcm,tgl,tgl)
			o.glcm[[k]] <- m.glcm/sum(m.glcm, na.rm=TRUE)
}
          o.glcm[[5]]<- (o.glcm[[1]]+o.glcm[[2]]+o.glcm[[3]]+o.glcm[[4]])/4
          names(o.glcm) <- c("th_0", "th_45", "th_90", "th_135", "ave")
return(list(glcm=o.glcm, level=t.level, d=d))
}}
########################################### 
# power spectram from fft() data  /swapping quadrants #
############################################
swap.quad <- function(x, disp=FALSE, reverse=FALSE) {
	p <- nrow(x); q <- ncol(x)
	if (reverse != TRUE) {
		if (p/2 == 0) {
		pp <- p %/% 2
	} else {
		pp <- ceiling(p/2)
	}
	if (q/2 == 0) {
		qq <- q %/% 2
	} else {
		qq <- ceiling(q/2)
	}
	} else {
		pp <- p %/% 2
		qq <- q %/% 2
	}
	ps <- matrix(0,p,q)
	x <- rbind(x,x)
	x <- cbind(x,x)
	ps <- x[1:p+pp, 1:q+qq]	
	if (disp == TRUE ){
	image(log(Mod(ps)), col=gray(c(0:255)/255), useRaster=T)		
	}
	return(ps)
}
#############################################
# polar transform from cartesian coodinates #
#############################################
car2pol <- function(x, method="bilinear") {
if (length(dim(x))>2){warning("data must be grayscale image")
} else{
if(max(x)<=1){
	x <- x*(2^attr(x,"bits.per.sample")-1)
}
	i.img <- x
	if (method =="NN") {
		p <-nrow(i.img)
		q <-ncol(i.img)
		m <- min(p,q)
		conv.img <- matrix(0,(trunc(m/2)-1),361)
		for ( th in 1:360) {# center is round(p/2), round(q/2)
			for (r in 1:(trunc(m/2)-1)) {
			th.rad <- pi*(th-1)/180
				conv.img[r,th] <- i.img[round(p/2+(r-1)*cos(th.rad))-1,round(q/2+(r-1)*sin(th.rad))]
			}
	}
	attr(conv.img, "bits.per.sample") <- attr(x,"bits.per.sample")
	attr(conv.img, "samples.per.pixel") <- attr(x, "samples.per.pixel")
	return(conv.img)
	} else if (method=="bilinear") {
		p <-nrow(i.img)
		q <-ncol(i.img)
		m <- min(p,q)
		xc <- p%/%2 ; yc <- q%/%2   # center of image
		conv.img <- matrix(0,(trunc(m/2)-1),361)
		for (th in 0:360+1) {# center is round(p/2), round(q/2)
			for (r in 1:(trunc(m/2)-1)) {
				th.rad <- pi*(th-1)/180
				Sxy <- c(xc+(r-1)*cos(th.rad), yc+(r-1)*sin(th.rad))
				xy <- trunc(Sxy)
				dxy <- Sxy-xy		
				a <- matrix(c(i.img[xy[1],xy[2]], i.img[xy[1]+1,xy[2]], i.img[xy[1],xy[2]+1],i.img[xy[1]+1,xy[2]+1]),2,2) 
				b <- matrix(c(1-dxy[1],dxy[1]),1,2) 
				c <- matrix(c(1-dxy[2],dxy[2]),2,1)
				conv.img[r,th] <- b %*% a %*% c
			}
		}}
	attr(conv.img, "bits.per.sample") <- attr(x,"bits.per.sample")
	attr(conv.img, "samples.per.pixel") <- attr(x, "samples.per.pixel")
	return(round(conv.img))
	if(method !="NN" && method!="bilinear"){
		return(NULL)
	}	
}}
##############################################
# Clockwise rotation of raster or matrix:    #
# Default is bilinear interpolationmethod    #
##############################################
rotate.matrix <- function(x, angle=10, method="bilinear"){
if (length(dim(x))>2){warning("data must be grayscale image")
} else{
if(max(x)<=1){
	x <- x*(2^attr(x,"bits.per.sample")-1)
}
	img <- x
	angle.rad <-angle*pi/180
	co.x <- matrix(rep(-(ncol(img)/2-0.5):(ncol(img)/2-0.5),nrow(img)),nrow=nrow(img),byrow=T)
	co.y <- matrix(rep(-(nrow(img)/2-0.5):(nrow(img)/2-0.5),ncol(img)),ncol=ncol(img))
	if(method=="simple"){
		co.xn <- round(co.x*cos(angle.rad)-co.y*sin(angle.rad)) # new coordinate
	co.yn <- round(co.x*sin(angle.rad)+co.y*cos(angle.rad))
	co.xn2 <- co.xn+max(co.xn)+1	# convert to topleft
	co.yn2 <- co.yn+max(co.yn)+1
	img.rot <-numeric(max(co.yn2)*max(co.xn2))
	img.rot[(co.xn2-1)*max(co.yn2)+co.yn2]<-img
	} else if(method=="NN" || method=="bilinear"){
	if(ncol(img)%%2==0){co.xn <- trunc(co.x*cos(angle.rad)-co.y*sin(angle.rad)+0.5) 
		}else{co.xn <- trunc(co.x*cos(angle.rad)-co.y*sin(angle.rad)) }
	if(nrow(img)%%2==0){co.yn <- trunc(co.x*sin(angle.rad)+co.y*cos(angle.rad)+0.5)
		}else{co.yn <- trunc(co.x*sin(angle.rad)+co.y*cos(angle.rad)) }
		# new coordinate
	co.list1 <- c(0,0)
	for(i in 1:max(co.xn)){
		y1 <- co.yn[which(co.xn==i)]
		y2 <- min(y1):max(y1)
		y3 <- matrix(c(rep(i,length(y2)),y2),nrow=2,byrow=T)
		co.list1 <- cbind(co.list1,y3)
	}
	co.list1 <- co.list1[,-1]
	co.list2 <-matrix(0,nrow(co.list1),ncol(co.list1))
	co.list2[1,] <- -co.list1[1,]
	if(ncol(img)%%2==0){co.list2[1,]<- co.list2[1,]+1}
	co.list2[2,] <- -co.list1[2,]
	if(nrow(img)%%2==0){co.list2[2,]<- co.list2[2,]+1}
	co.list <- cbind(co.list1,co.list2)
	if(ncol(img)%%2!=0){
		y0 <- co.yn[which(co.xn==0)]
		y02 <- min(y0):max(y0)
		y03 <- matrix(c(rep(0,length(y02)),y02),nrow=2,byrow=T)
		co.list <- cbind(co.list,y03)
	}
	co.xn2 <- co.list[1,]+max(co.list[1,])+1	# convert to topleft
	co.yn2 <- co.list[2,]+max(co.list[2,])+1
	if(ncol(img)%%2==0){co.list[1,] <- co.list[1,] -0.5}
	if(nrow(img)%%2==0){co.list[2,] <- co.list[2,] -0.5}
	if(method=="NN"){
		if(ncol(img)%%2==0){co.xn.b <- round(co.list[1,]*cos(-angle.rad)-co.list[2,]*sin(-angle.rad)+0.5)+ncol(img)/2
			}else{co.xn.b <- round(co.list[1,]*cos(-angle.rad)-co.list[2,]*sin(-angle.rad))+ncol(img)/2+0.5}
		if(nrow(img)%%2==0){co.yn.b <- round(co.list[1,]*sin(-angle.rad)+co.list[2,]*cos(-angle.rad)+0.5)+nrow(img)/2
			}else{co.yn.b <- round(co.list[1,]*sin(-angle.rad)+co.list[2,]*cos(-angle.rad))+nrow(img)/2+0.5}
		img.rot <-numeric(max(co.yn2)*max(co.xn2))
		img.rot[(co.xn2-1)*max(co.yn2)+co.yn2]<-img[(co.xn.b-1)*max(co.yn.b)+co.yn.b]
	} else if(method=="bilinear"){
		co.xn.b <- co.list[1,]*cos(-angle.rad)-co.list[2,]*sin(-angle.rad)+0.5+ncol(img)/2
		co.yn.b <- co.list[1,]*sin(-angle.rad)+co.list[2,]*cos(-angle.rad)+0.5+nrow(img)/2
		co.x0 <- floor(co.xn.b)+1
		co.y0 <- floor(co.yn.b)+1
		ax <- co.xn.b-co.x0+1
		ay <- co.yn.b-co.y0+1
		img.e<-cbind(img[,1],img,img[,ncol(img)])
		img.e <-rbind(img.e[1,],img.e,img.e[nrow(img.e),])
		img.rot <-numeric(max(co.yn2)*max(co.xn2))
		img.rot[(co.xn2-1)*max(co.yn2)+co.yn2]<-(1-ax)*(1-ay)*img.e[(co.x0-1)*nrow(img.e)+co.y0]+ax*(1-ay)*img.e[co.x0*nrow(img.e)+co.y0]+ay*(1-ax)*img.e[(co.x0-1)*nrow(img.e)+co.y0+1]+ax*ay*img.e[co.x0*nrow(img.e)+co.y0+1]
	}
	}
	dim(img.rot)<- c(max(co.yn2),max(co.xn2))
	attr(img.rot, "bits.per.sample") <- attr(img,"bits.per.sample")
	attr(img.rot, "samples.per.pixel") <- attr(img, "samples.per.pixel")
	return(round(img.rot))
}}
###############################
# Haralick texture parameters              #
###############################
haralick <- function(x) {
o.hara <- matrix(0, 15,5) #Haralick parameters output matrix
for (th in 1:5) {
pglcm<-x$glcm[[th]]
nx <- ncol(pglcm)
ny <- nrow(pglcm)
px <- colSums(pglcm)
py <- rowSums(pglcm)
pxpy <-matrix(px,nx,ny)*t(matrix(py,nx,ny))
px_y <- matrix (0, nx+ny)
pxmy <- matrix (0, (nx+ny)/2)

# means & standard deviation
vx <- 1:nx
vy <- 1:ny
mx <- sum(px*vx) 
my <- sum(py*vx)
stdevx <- sum(px*(vx-mx)^2)
stdevy <- sum(py*(vy-my)^2)

# HX,HY,HXY for f12 and f13
hxy1_0 <- matrix (0, nx,ny)
hxy2_0 <- matrix (0, nx,ny)
hxy1_0 <- pglcm*log10(pxpy)
hxy2_0 <- (pxpy)*log10(pxpy)

hx <- -sum(px*log10(px),na.rm=TRUE)
hy <- -sum(py*log10(py),na.rm=TRUE)
hxy1 <- -sum(hxy1_0,na.rm=TRUE)
hxy2 <- -sum(hxy2_0,na.rm=TRUE)
op <- matrix(1:nx,nx,ny)
oq <- t(op)
spq <- matrix(1:nx,nx,ny)+t(matrix(1:ny,nx,ny))
dpq <- abs(matrix(1:nx,nx,ny)-t(matrix(1:ny,nx,ny)))

#1 Angular Second Moment / Homogeniety "asm"
o.hara[1,th] <- sum(pglcm^2)				
#2 Contrast "con"
o.hara[2,th] <- sum(dpq^2*pglcm)
#3 inverse Difference Moment "idm"
o.hara[3,th] <- sum(pglcm/(1+dpq^2))
#4 Entropy "ent"
o.hara[4,th] <- -sum(pglcm*log10(pglcm),na.rm=TRUE) 
#5 Correlation 	"cor"
o.hara[5,th] <- sum ((op-mx)*(oq-my)*pglcm/(sqrt(stdevx*stdevy)))
#6 Variance in Haralick 1973	"var"
o.hara[6,th] <- sum((op-((mx+my)/2))^2*pglcm)
#7 Sum Average "sav"
o.hara[7,th] <- sum(spq*pglcm) 
#8 Sum Entropy "sen"
#9 Difference Entropy "den"
sen<- array(0,(2*nx))  # sen
den.1 <- array(0,nx)  # den
den.2 <- array(0,nx)  # den
pglcm2 <- cbind(pglcm[,nx:1])  # a matrix with its column reverse order
for (i in 2:nx) {
sen[i]<-sum(diag(pglcm2[1:i,(nx-i+1):nx]))  # sen upper diagonal (include diagonal)
den.1[i]<-sum(diag(pglcm[1:i,(nx-i+1):nx]))  # den lower diagonal (include diagonal)
sen[1]<-pglcm2[1,nx]
den.1[1]<-pglcm[1,nx]
}
for (i in 1:(nx-2)) {
sen[i+nx]<-sum(diag(pglcm2[(i+1):nx,1:(nx-i)]))  # sen upper diagonal (include diagonal)
den.2[nx-i]<-sum(diag(pglcm[(i+1):nx,1:(nx-i)]))  # den lower diagonal (include diagonal)
}
sen[nx+nx-1]<-pglcm2[nx,1]
den.2[1]<-pglcm[nx,1]
o.hara[8,th] <- -sum(sen*log10(sen),na.rm=TRUE)  # sen
den <- den.1+den.2
o.hara[9,th] <- -sum(den*log10(den),na.rm=TRUE)  # den
#10 Difference Variance "dva"
o.hara[10,th]<- sum(((dpq-o.hara[9])^2)*pglcm )
#11 Sum Variance "sva"
o.hara[11,th] <- sum(((spq-o.hara[8])^2)*pglcm)
#12 Information Measures of Correlation "f12" (- sign was intentionally added as the value give minus value)
o.hara[12,th] <- -(o.hara[4]-hxy1)/max(hx,hy)
#13 Information Measures of Correlation "f13"
o.hara[13,th] <- sqrt(1-exp(-2*abs(hxy2-o.hara[4] ))) 
#14 Cluster Shade "sha"
o.hara[14,th] <- sum((spq-mx-my)^3*pglcm)
#15 Cluster prominence "pro"
o.hara[15,th] <- sum((spq-mx-my)^4*pglcm)	#15 Cluster prominence
	}
	colnames(o.hara) <- c("th_0", "th_45", "th_90", "th_135", "ave")
	rownames(o.hara) <- c("asm","con","idm", "ent","cor","var","sav","sen","den","dva","sva","f12","f13", "sha","pro")
	rng <- apply(o.hara,1,max)-apply(o.hara,1,min)
     o.hara <- cbind(o.hara,rng)
return(o.hara)
} #  (15 parameters calculated)

#########################################
###matrix rotation 90 degree clockwise###
#########################################
rot90c <- function(x){
	img.rot <- t(apply(x,2,rev))
	attr(img.rot, "bits.per.sample") <- attr(x,"bits.per.sample")
	attr(img.rot, "samples.per.pixel") <- attr(x, "samples.per.pixel")
	return(img.rot)
	
}
######################################
###Gabor Filter in Frequency Domain###
######################################
gabor.filter <- function(x, lamda=5, theta=45, bw=1.5, phi=0, asp=1, disp=FALSE ) {
if (length(dim(x))>2){warning("data must be grayscale image")
} else {
img <- x
if (min(range(dim(img))) > 151){ 		#Gabor Kernel 
	ks <- 151
} else {
	ks <- min(range(dim(img)))-5
}	
	rot <- pi*theta/180						#rotation angle in radian       
	phi <- pi*phi/180						#offset phase angle in radian
	sigma <- sqrt(log(2,2)/2)*(2^bw+1)*lamda/(pi*(2^bw-1))	#band width to sdev of envelope
if (ks %% 2 == 0){
	ks <- ks+1
	m <-trunc(ks %/% 2)	
} else {
	m <- trunc(ks %/% 2)
}
	gf_Re <- matrix(0,ks,ks)
	gf_Im <- matrix(0,ks,ks)
	for ( x in -m:m ) {
		for ( y in -m:m ) {
			rx <- x*cos(rot)+y*sin(rot)  											# rotation of  x coodinate
			ry <- -x*sin(rot)+y*cos(rot) 											# rotation of  y coodinate
			gauss <- exp((-rx^2-(ry*asp)^2)/(2*sigma^2))				# envelope function
			gf_Re[m+x+1,m+y+1] <- gauss*cos(2*rx*pi/lamda+phi)			# Gabor filter in real number
			gf_Im[m+x+1,m+y+1] <- gauss*sin(2*rx*pi/lamda+phi)				# Gabor filter in imaginary number
		}
	}
kernel <- matrix(complex(real=gf_Re, imaginary=gf_Im), ks, ks)
#
mx <- nrow(img) ; my <- ncol(img)  			# size of an imput image to be filtered
cx <- trunc(mx/2)									# x center of gabor filter
cy <- trunc(my/2)									# y center of gabor filter
mask <- matrix(0+0i,mx,my)					# filter mask for inverse FFT
mask[-m:m+cx,-m:m+cy] <- gf_Re			# locate Gabor fuction in the center
fd_mask <- fft(Re(mask))						# transform of the real part of Gabor filter in frequency domain
img_out <- fft(Mod(fd_mask)*fft(img), inverse=TRUE)  # convolution in frequency domain and inverse FFT
#
if (disp == TRUE) {
	par(mfrow=c(2,2),mar=c(3,3,3,3))
	comment <- paste("lmd", lamda, "th", theta, "ph", phi, "bw", bw, "asp", asp, sep="")
	image(rot90c(img), asp=1, main="original", col=gray(c(0:255)/255), axes=FALSE, useRaster=TRUE)
	image(rot90c(Re(mask)), asp=1, main=comment, col=gray(c(0:255)/255), axes=FALSE, useRaster=TRUE)
	image(rot90c(Mod(swap.quad(fd_mask))), asp=1, main="mask in frequency domain", col=gray(c(0:255)/255), axes=FALSE, useRaster=TRUE)
	image(rot90c(Mod(img_out)), asp=1, main="filtered", col=gray(c(0:255)/255), axes=FALSE, useRaster=TRUE)
	return(list(kernel=kernel, mask=mask, freq_mask=fd_mask, filtered_img=Mod(img_out)))
} else {
	return(list(kernel=kernel, mask=mask, freq_mask=fd_mask, filtered_img=Mod(img_out)))
}
}}
######################################
###Edge Detection
######################################
edge.detect <- function(x, thresh1=1, thresh2=15, noise="gaussian", noise.s=3, method="Canny") {
if (length(dim(x))>2){warning("data must be grayscale image")
} else {
img <- x
img.n <- noise.filter(img,noise.s, method= noise)
img.ed <- img[2:(nrow(img)-1),2:(ncol(img)-1)]
rasSX <- matrix(NA,length(img.ed),9)
rasSY <- matrix(NA,length(img.ed),9)
sx <- c(-1,-2,-1,0,0,0,1,2,1)/8
sy <- c(-1,0,1,-2,0,2,-1,0,1)/8
for(i in 1:3){
	for(j in 1:3){
		rasSX[,(3*(j-1)+i)] <- sx[3*(j-1)+i]*array(img.n[(i:(nrow(img.n)+i-3)),(j:(ncol(img.n)+j-3))])
		rasSY[,(3*(j-1)+i)] <- sy[3*(j-1)+i]*array(img.n[(i:(nrow(img.n)+i-3)),(j:(ncol(img.n)+j-3))])
		}
	}
img.sx <- apply(rasSX,1,sum)
img.sy <- apply(rasSY,1,sum)
dim(img.sx) <- dim(img.sy) <- dim(img.ed)
img.sxsy <- (img.sx^2+img.sy^2)^0.5			   #magnitude of gradient
if(method=="Sobel"){
	out <- img.sxsy
} else if(method=="Canny"){
	img.th <- atan(img.sy/img.sx)
	img.th0 <- img.th						       #angles of gradient (0,45,90,135)
	img.th0[abs(img.th0) >= pi*3/8] <- 90
	img.th0[abs(img.th0) <= pi/8] <- 0
	img.th0[img.th0 > pi/8 & img.th0 < pi*3/8] <- 45
	img.th0[img.th0 > -pi*3/8 & img.th0 < -pi/8] <- 135
	#non-maximum suppression
	img.sxsy.bl <- matrix(0,nrow(img),ncol(img))
	img.sxsy.bl[2:(nrow(img)-1),2:(ncol(img)-1)] <- img.sxsy
	rasL <- matrix(NA,length(img.sxsy),9)
	for(i in 1:3){
	for(j in 1:3){
		rasL[,(3*(j-1)+i)] <- array(img.sxsy.bl[(i:(nrow(img.n)+i-3)),(j:(ncol(img.n)+j-3))])
		}
	}
	rasLj <- matrix(0,length(img.sxsy),4)
	rasLj[which(rasL[,5] > rasL[,4] & rasL[,5] > rasL[,6]),3] <- 1
	rasLj[which(rasL[,5] > rasL[,7] & rasL[,5] > rasL[,3]),4] <- 1
	rasLj[which(rasL[,5] > rasL[,2] & rasL[,5] > rasL[,8]),1] <- 1
	rasLj[which(rasL[,5] > rasL[,1] & rasL[,5] > rasL[,9]),2] <- 1
	img.tl <- array(img.sxsy)
	img.tl[which(img.th0==0)] <- img.sxsy[which(img.th0==0)]*rasLj[which(img.th0==0),1]
	img.tl[which(img.th0==45)] <- img.sxsy[which(img.th0==45)]*rasLj[which(img.th0==45),2]
	img.tl[which(img.th0==90)] <- img.sxsy[which(img.th0==90)]*rasLj[which(img.th0==90),3]
	img.tl[which(img.th0==135)] <- img.sxsy[which(img.th0==135)]*rasLj[which(img.th0==135),4]
	dim(img.tl) <- dim(img.sxsy)
	#hysteresis threshold
	mxth <- sort(img.tl[which(img.tl!=0)])[round(length(img.tl[which(img.tl!=0)])*thresh2/100)]
	mnth <- sort(img.tl[which(img.tl!=0)])[round(length(img.tl[which(img.tl!=0)])*thresh1/100)]
	img.bn0 <- img.bn1 <- matrix(0,nrow(img.tl)+2,ncol(img.tl)+2)
	img.bn0[2:(nrow(img.bn0)-1),2:(ncol(img.bn0)-1)] <- img.bn1[2:(nrow(img.bn0)-1),2:(ncol(img.bn0)-1)] <- img.tl
	img.bn0[which(img.bn0 <= mnth)] <- 0
	img.bn0[which(img.bn0 > mnth)] <- 1	
	img.lb <- matrix(0,nrow(img.tl)+2,ncol(img.tl)+2)
	ren8 <- rbind(c(-1,-1),c(-1,0),c(-1,1),c(0,-1)) #8-connected area
	k <- 0
	for(i in 2:(nrow(img.bn0)-1)){
	for(j in 2:(ncol(img.bn0)-1)){
		z8 <- t(c(i,j)+t(ren8))
		if(img.bn0[i,j]==0){
			}else if(img.bn0[i,j]==1 && sum(img.bn0[z8]==0)==4){
				k <- k+1
				img.lb[i,j] <- k
			}else {
				zl <- z8[which (img.lb[z8] != 0),,drop=F]
				z.lb <- img.lb[zl]
				z.lb <- sort(z.lb[!duplicated(z.lb)])
				img.lb[i,j] <- z.lb[1]
				if(length(z.lb)>1){
					for(l in 2:length(z.lb)){
				 	img.lb[which(img.lb==z.lb[l])] <- z.lb[1]
						}
					}

			}
		}
		}
	img.bn <- img.lb
	lb <- sort(img.lb[!duplicated(array(img.lb))])
	lb <- lb[-1]
	for(i in lb){
		if(length(which(img.bn1[which(img.lb==i)]>mxth))==0){
			img.bn[which(img.lb==i)] <- 0
		}else {
			img.bn[which(img.lb==i)] <- 1
	}
	}
	img.bn <- img.bn[2:(nrow(img.bn0)-1),2:(ncol(img.bn0)-1)]
	out <- img.bn
}else {
	print(NULL)
}
}}
######################################
###Connected Component Labelling
######################################
cc.label <- function(x, connect=8, inv=FALSE, img.show=FALSE,text.size=0.3){
if (length(x[which(x==0 | x==1)])!=length(x)){
	print("x should be a binary image.")
} else {
	if(inv){
		x <- abs(x-1)
	}
	img <- img.lb <- matrix(0,nrow(x)+2,ncol(x)+2)
	img[2:(nrow(img)-1),2:(ncol(img)-1)] <- x
	if(connect==8){
		conn <- rbind(c(-1,-1),c(-1,0),c(-1,1),c(0,-1)) #8-connected area
	} else if(connect==4){
		conn <- rbind(c(-1,0),c(0,-1)) #4-connected area
	} else {
		print("connect should be 4 or 8")
	}	
	
	k <- 0
	for(i in 2:(nrow(img.lb)-1)){
	for(j in 2:(ncol(img.lb)-1)){
		z <- t(c(i,j)+t(conn))
		if(img[i,j]==0){
			}else if(img[i,j]==1 && sum(img[z]==0)==nrow(conn)){
				k <- k+1
				img.lb[i,j] <- k
			}else {
				zl <- z[which (img.lb[z] != 0),,drop=F]
				z.lb <- img.lb[zl]
				z.lb <- sort(z.lb[!duplicated(z.lb)])
				img.lb[i,j] <- z.lb[1]
			if(length(z.lb)>1){
				for(l in 2:length(z.lb)){
				 img.lb[which(img.lb==z.lb[l])] <- z.lb[1]
					}
				}
			}
		}
		}
lb <- sort(img.lb[!duplicated(array(img.lb))])
lb <- lb[-1]
for(i in lb){
		img.lb[which(img.lb==i)] <- which(lb==i)
	}
img.lb <- img.lb[2:(nrow(img.lb)-1),2:(ncol(img.lb)-1)]
n.lb <- length(lb)
sum.lab <- data.frame(matrix(0,n.lb,7))
colnames(sum.lab) <- c("label","area","aveX","aveY","dX","dY","edge")
sum.lab[,1] <- 1:n.lb
for (i in 1:n.lb){
	cdn <- which(img.lb==i,arr.ind=T)
	if (length(cdn[which(cdn==1 | cdn[,1]== nrow(img.lb) | cdn[,2]==ncol(img.lb))]) ==0){edg <- 0
		}else {edg <- 1}
	ranX <- range(cdn[,1])
	ranY <- range(cdn[,2])
	sum.lab[i,2:7] <- c(nrow(cdn),mean(cdn[,1]),mean(cdn[,2]),ranX[2]-ranX[1]+1,ranY[2]-ranY[1]+1,edg)
	}
if(img.show){
		image(rot90c(img),col=(255:0)/255,ann=F,axes=F,useRaster=TRUE)
		par(new=T)
		plot(sum.lab$aveY,sum.lab$aveX,type="n",ylim=c(nrow(img),1),xlim=c(1,ncol(img)),xaxs="i",yaxs="i",ann=F,axes=F)
		text(sum.lab$aveY,sum.lab$aveX,sum.lab$label,cex=text.size,ylim=c(nrow(img),1),xlim=c(1,ncol(img)),col="red",xaxs="i",yaxs="i",ann=F)
	}
	return(list(image=img.lb,summary=sum.lab))
}
}
######################################
###Cropping
######################################
crop <- function(x, width=300, height=300, shift=c(0,0)){
if (length(dim(x))>2){warning("data must be grayscale image")
} else {img <- x
	c.x <- nrow(img)/2+0.5
	c.y <- ncol(img)/2+0.5
	tl.x <- trunc(c.x-height/2+0.5)+shift[1]
	tl.y <- trunc(c.x-width/2+0.5)+shift[2]
	img.o <- img[tl.x:(tl.x+height-1),tl.y:(tl.y+width-1)] 
	attr(img.o, "bits.per.sample") <- attr(img,"bits.per.sample")
	attr(img.o, "samples.per.pixel") <- attr(img, "samples.per.pixel")
	return(img.o)
}
}