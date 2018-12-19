################################################################################
# PAHO Mortality Data                                                          #
#                                                                              #
#        DATE: October 2018                                                    #
#    CODED BY: Kayoko Shioda (kayoko.shioda@yale.edu)                          #
#     ADVISOR: Dan Weinberger (daniel.weinberger@yale.edu)                     #
################################################################################

#------------------------------------------------------------------------------#
# DESCRIPTION
#------------------------------------------------------------------------------#

# Figure (Forest plot)
# RR estimated by the "best" model, stratified by age group and country

#------------------------------------------------------------------------------#
# Load datasets
#------------------------------------------------------------------------------#

rm(list = ls(all = TRUE))


#select what you want to plot : "best" estimate; full (sc) estimates, time-trend, or stl_pca
ds.select<-'rr_best'  #rr_best, rr_full, rr_time_trend, rr_pca


# Load dataset for each country and merge 10 files to create one big file.
#ar_SubChpt_BestModel_woPandemic_rr_full.csv
#ar_SubChpt_BestModel_woPandemic_rr_full.csv
################################################################################
# For HDI-level estimates without unbiasing
setwd("C:/Users/dmw63/Weinberger Lab Dropbox/PAHO mortality/Results/HDI/SubChapterLevel/Monthly/BEST MODEL code results wo Pandemic")
country <- "ar"
dt <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt) <- c("age_group","Lower","Median","Upper")
country <- "br"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "co"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "dr"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "ec"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "gy"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "hr"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "mx"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_rr_best_National.csv",sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "nc"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)
country <- "pr"
dt1 <-read.csv(file = paste(country,"_SubChpt_BestModel_woPandemic_", ds.select, '.csv',sep=""))
colnames(dt1) <- c("age_group","Lower","Median","Upper")
dt <- rbind(dt,dt1)

# ################################################################################
# # For HDI-level estimates **with** unbiasing
# setwd("/Users/shiodakayoko/Weinberger Lab Dropbox/PAHO mortality/Results/HDI/SubChapterLevel/Monthly/best model wo pandemic unbias")
# country <- "PAHO_ar"
# dt <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt) <- c("age_group","Lower","Median","Upper")
# country <- "PAHO_br"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_co"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_dr"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_ec"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_gy"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_hr"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_mx"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_nc"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)
# country <- "PAHO_pr"
# dt1 <-read.csv(file = paste(country,"/",country,"_", ds.select, '.csv',sep=""))
# colnames(dt1) <- c("age_group","Lower","Median","Upper")
# dt <- rbind(dt,dt1)

#------------------------------------------------------------------------------#
# Set up
#------------------------------------------------------------------------------#

dt$eb_l <- dt$Median - dt$Lower
dt$eb_u <- dt$Upper - dt$Median

# Let's exclude HDI-level estimates now...
dt <- dt[substr(dt$age_group,
                nchar(as.character(dt$age_group)), 
                nchar(as.character(dt$age_group)))=="A",]

# Create "age" and "agegrp"
dt$age <- NA
for (i in 1:nrow(dt)) {
  dt$age[i] <- as.character(unlist(strsplit(as.character(dt$age_group[i]), " ", fixed = TRUE))[2])
}
table(dt$age)
dt$agegrp <- ifelse(dt$age=="<2m", 0, NA)
dt$agegrp <- ifelse(dt$age=="2-11m",  1, dt$agegrp)
dt$agegrp <- ifelse(dt$age=="12-23m", 2, dt$agegrp)
dt$agegrp <- ifelse(dt$age=="24-59m", 3, dt$agegrp)
dt$agegrp <- ifelse(dt$age=="2-23m",  4, dt$agegrp)
dt$agegrp <- ifelse(dt$age=="2-59m",  5, dt$agegrp)
table(dt$age,dt$agegrp)

# Create "country"
dt$country <- NA
for (i in 1:nrow(dt)) {
  dt$country[i] <- as.character(unlist(strsplit(as.character(dt$age_group[i]), " ", fixed = TRUE))[1])
}


#------------------------------------------------------------------------------#
# Plots
#------------------------------------------------------------------------------#

# #-----*-----*-----*-----*-----*-----*-----#
# # RR for one selected stratum
# #-----*-----*-----*-----*-----*-----*-----#
# 
# # Create a subset of the data
# sub <- dt[which(dt$age=="2-59m"),]
# 
# # Graph
# library(ggplot2)
# graph <- ggplot(sub, aes(x=Median, y=country)) + 
#   geom_errorbarh(aes(xmin=Median-eb_l, xmax=Median+eb_u), colour="#969696", height=0.1) +
#   geom_line(linetype="blank") +geom_point(size=1) +
#   scale_x_continuous(name="",limits=c(0,2)) + xlab("") + ylab("") +
#   theme(text = element_text(size=20),panel.background = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         axis.line.x = element_line(color="black", size = 0.3),
#         axis.line.y = element_blank(),axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()) +
#   geom_vline(xintercept=c(1),lty=c(2),color=c("darkgrey"))
# graph
# #zoom <- graph + coord_cartesian(xlim=c(0,2.5))


#-----*-----*-----*-----*-----*-----*-----#
# RR by age group and country
#-----*-----*-----*-----*-----*-----*-----#

# Sort the dataset by age group and country
dt <- dt[order(dt$agegrp,dt$country),]

# Insert empty rows to have some space between each age group
newrow <- matrix(rep(NA,4*ncol(dt)), nrow=4, ncol=ncol(dt))
colnames(newrow) <- colnames(dt)

dt<-dt[dt$age=='2-59m' | dt$age =='2-23m',] #Exclude the aggregated groups
dt$country.grp<-NA
dt$country.grp[dt$country %in% c('ar','br','co','pr','mx')]<-3
dt$country.grp[dt$country %in% c('ec','hr','nc')]<-2
dt$country.grp[dt$country %in% c('dr','gy')]<-1
#Country labels
dt$country.name<-NA
dt$country.name[dt$country=='ar']<-'Argentina'
dt$country.name[dt$country=='br']<-'Brazil'
dt$country.name[dt$country=='co']<-'Colombia'
dt$country.name[dt$country=='dr']<-'Dom. Rep.'
dt$country.name[dt$country=='ec']<-'Ecuador'
dt$country.name[dt$country=='hr']<-'Honduras'
dt$country.name[dt$country=='nc']<-'Nicaragua'
dt$country.name[dt$country=='pr']<-'Peru'
dt$country.name[dt$country=='mx']<-'Mexico'
dt$country.name[dt$country=='gy']<-'Guyana'

#Metafor analysis and plotting
#install.packages('metafor')
library(metafor)
age.list<-c( '2-23m','2-59m')
par(mfrow=c(1,1), mar=c(4,1,1,1))
for(i in 1: length(age.list)){
age.plot<-dt[dt$age==age.list[i],]
age.plot<-age.plot[!is.na(age.plot$Median),]
age.plot$log.rr<-log(age.plot$Median)
age.plot$approx.var<- ((log(age.plot$Upper) - log(age.plot$Median))/1.96)^2
ma<-rma.uni( yi=log.rr, vi=approx.var, data=age.plot, slab=country.name )
forest(ma,atransf=exp)
}

#http://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups
overall.plot<-dt[!is.na(dt$Median),]
# ma.test<-rma.uni( yi=log.rr, vi=approx.var,mods = ~  age ,data=overall.plot, slab=country.name )
# forest(ma.test,atransf=exp)

#ma.all<-rma.uni( yi=log.rr, vi=approx.var, data=overall.plot, slab=country )
overall.plot$log.rr<-log(overall.plot$Median)
overall.plot$approx.var<- ((log(overall.plot$Upper) - log(overall.plot$Median))/1.96)^2
ma.2_23m<-rma.uni( yi=log.rr, vi=approx.var, data=overall.plot, slab=country , subset=(age=='2-23m'))
ma.2_59m<-rma.uni( yi=log.rr, vi=approx.var, data=overall.plot, slab=country , subset=(age=='2-59m'))

#simple plot, unstratified
#par(mfrow=c(1,1), mar=c(1,1,1,1))
#forest(x=overall.plot$Median, ci.lb=overall.plot$Lower, refline=1, ci.ub=overall.plot$Upper, slab=overall.plot$country)

#Stratified
#Split by age group to calculate number of data points for each
plot.spl<-split(overall.plot, overall.plot$agegrp)
plot.spl<-lapply(plot.spl, function(x) x<-x[order(x$country.grp, x$country),])
n.country.grp<-sapply(plot.spl, function(x) nrow(x) )
#Specify where each group starts and ends, providing space between groups
start.grp1<-3
start.grp2<-start.grp1+n.country.grp[2]+4

plot.indices<-c(start.grp1:(start.grp1+n.country.grp[2]-1) ,
                start.grp2:(start.grp2+n.country.grp[1]-1)
                )
overall.plot.sorted<-overall.plot[order(-overall.plot$agegrp, -overall.plot$country.grp, overall.plot$country),]
### decrease margins so the full space is used
tiff(paste0('C:/Users/dmw63/Weinberger Lab Dropbox/PAHO mortality/Results/Forest plots/agg without unbiasing/forest.natl_agg',ds.select,'.tiff'),
        width=5, height=6, units='in', res=200)
par(mar=c(4,4,1,2))
forest(x=overall.plot.sorted$Median,rows=plot.indices, 
       ci.lb=overall.plot.sorted$Lower, 
       vi=overall.plot.sorted$approx.var,
        ci.ub=overall.plot.sorted$Upper,
       refline=1, 
          slab=overall.plot.sorted$country.name,
            ylim=c(0,(max(plot.indices)+4.5)), clim=c(0.4, 1.5),xlim=c(0, 2), 
          at=c(0.4,0.6,0.8,1,1.2,1.4,1.6),
            cex=0.75)
### add summary polygons for the three subgroups
#addpoly(ma.2_59m, row=start.grp1-1.5, cex=0.75, transf=exp, mlab="")
#addpoly(ma.2_23m, row= start.grp2-1.5, cex=0.75, transf=exp, mlab="")

### add text for the subgroups
### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.75, font=4)
text(0, c((start.grp2+n.country.grp[1]+0.5),
            (start.grp1+n.country.grp[2]+0.5)), 
            pos=4, c("2-23m",'2-59m'), cex=0.9)
### switch to bold fonts
par(font=2)
text(0,                (max(plot.indices)+4), "Country",  pos=4)
text(2, (max(plot.indices)+4), "Rate Ratio [95% CrI]", pos=2)

#text(0,start.grp1-1.5, 'RE Model 2-59m' ,  pos=4,cex=0.75)
#text(0,start.grp2-1.5, 'RE Model 2-23m' ,  pos=4, cex=0.75)
dev.off()


#Playingaroun : is there a relationship between completeness of death registry and reported decline (expect bigger declines with less complete registries)
#https://data.worldbank.org/indicator/SP.REG.DTHS.ZS?end=2012&locations=BR-DO-AR-EC-GY-HN-MX-NI-CO-PE&name_desc=false&start=1990&view=chart&year=2007
overall.plot$completeness<-NA
overall.plot$completeness[overall.plot$country=='ar']<-100
overall.plot$completeness[overall.plot$country=='br']<-93
overall.plot$completeness[overall.plot$country=='co']<-98
overall.plot$completeness[overall.plot$country=='dr']<-57
overall.plot$completeness[overall.plot$country=='ec']<-80
overall.plot$completeness[overall.plot$country=='hr']<-17
overall.plot$completeness[overall.plot$country=='nc']<-68
overall.plot$completeness[overall.plot$country=='mx']<-99
overall.plot$completeness[overall.plot$country=='pe']<-69
overall.plot$completeness[overall.plot$country=='gy']<- 81
ages<-unique(overall.plot$age)
par(mfrow=c(1,2))
for(i in ages){
plot(overall.plot$completeness[overall.plot$age==i], overall.plot$log.rr[overall.plot$age==i], bty='l',ylim=c(-0.5,0.5))
  abline(h=0)
  title(i)
}