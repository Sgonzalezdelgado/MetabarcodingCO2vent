
######## Metabarcoding CO2vent analysis ########

  # Authors: Owen S. Wangensteen, owenwangensteen@ub.edu; Sara González-Delgado, sgonzadl@ull.edu.es
  
library(vegan)
library(readxl)
library(readr)
library(dunn.test)

## Funcion renormalize
renormalize <- function(db){
  total_reads <- colSums(db)
  for (i in 1:ncol(db)) db[,i] <- db[,i]/total_reads[i]
  return(db)
}

# Load the databases. 
  #Remember to create columns with the groupings of interest (e.g= group1 (Algae, Metzoans), group3 (Calcifiers, Non-calcifiers).

db_swarm <- read_delim("~/Desktop/metabarcodinganalysis/VENT-SWARM13-algaemetazoan.csv", 
                                        delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Check the column names
colnames(db_swarm)

# Define function for calculating relative abundances per sample
normalize_samples <- function(db,sample_cols){
  sumcols <- colSums(db[,sample_cols])
  for (i in sample_cols) db[,i] <- db[,i]/sumcols[i+1-sample_cols[1]]
  return(db)
}

# Calculate normalized databases
db_swarm_normalized <- normalize_samples(db_swarm,17:64)
View(db_swarm_normalized)

#### Alpha diversity ####

# Rarefaction curve
db <- db_swarm
sample_columns <- 17:64

#for (i in sample_columns) db[,i]<- as.numeric(db[,i])
abundances <- db[,sample_columns]
total_abundances <- rowSums(abundances)
abundances <- abundances[total_abundances>1,]
names(abundances)
total_reads <- colSums(abundances)
max(total_reads)
region <- paste0(substr(colnames(abundances),2,2),substr(colnames(abundances),4,4))
region <- as.factor(region)

colors=c("red3","orange","palegreen","cornflowerblue", "red4", "tomato2", "olivedrab", "darkblue")

plot_rarecurve <- function(samp_abun,main="Rarefaction curves",xmax=60000,ymax=250){
  rarecurve(t(samp_abun),step=50,label=F,col=colors[1:ncol(samp_abun)],lwd=1,las=1,xlim=c(0,xmax),ylim=c(0,ymax),main=main,xlab="Number of reads",ylab="Number of MOTUS")
  legend("bottomright",names(samp_abun),lty = 1,col=colors[1:ncol(samp_abun)],cex=0.5,bty = "n")
}

pdf("Rarefaction plot all Vent COI Swarm.pdf", width = 8, height = 11)
par(mfrow=c(3,2),oma=c(2,2,2,2))
plot_rarecurve(abundances[,region=="1A"],"Vent A")
plot_rarecurve(abundances[,region=="1B"],"Vent B")
plot_rarecurve(abundances[,region=="2A"],"Transition25 A")
plot_rarecurve(abundances[,region=="2B"],"Transition25 B")
plot_rarecurve(abundances[,region=="3A"],"Transition75 A")
plot_rarecurve(abundances[,region=="3B"],"Transition75 B")
plot_rarecurve(abundances[,region=="4A"],"Control A")
plot_rarecurve(abundances[,region=="4B"],"Control B")
dev.off()

#Define function for plotting barplots and calculating Dunn's tests

# RICHNESS
rarefaction_dunn_sr <- function(db,name_output){
  abundances <- db[,(substr(names(db),1,1)=="V")]
  rownames(abundances) <- as.character(paste(db$scientific_name,db$id,sep="_"))
  motu_reads <-  rowSums(abundances)
  abundances <- abundances[order(motu_reads,decreasing = T),]
  abundances <- abundances[rowSums(abundances)>0,]
  # Species richness by rarefaction
  repl.raref <- replicate(50,apply(rrarefy(t(abundances),sample=5000)>0,1,sum),simplify="array")
  viol_V11A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="11A"),])
  viol_V12A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="12A"),])
  viol_V13A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="13A"),])
  viol_V14A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="14A"),])
  viol_V15A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="15A"),])
  viol_V16A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="16A"),])
  viol_V11B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="11B"),])
  viol_V12B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="12B"),])
  viol_V13B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="13B"),])
  viol_V14B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="14B"),])
  viol_V15B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="15B"),])
  viol_V16B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="16B"),])
  viol_V21A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="21A"),])
  viol_V22A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="22A"),])
  viol_V23A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="23A"),])
  viol_V24A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="24A"),])
  viol_V25A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="25A"),])
  viol_V26A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="26A"),])
  viol_V21B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="21B"),])
  viol_V22B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="22B"),])
  viol_V23B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="23B"),])
  viol_V24B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="24B"),])
  viol_V25B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="25B"),])
  viol_V26B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="26B"),])
  viol_V31A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="31A"),])
  viol_V32A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="32A"),])
  viol_V33A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="33A"),])
  viol_V34A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="34A"),])
  viol_V35A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="35A"),])
  viol_V36A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="36A"),])
  viol_V31B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="31B"),])
  viol_V32B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="32B"),])
  viol_V33B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="33B"),])
  viol_V34B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="34B"),])
  viol_V35B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="35B"),])
  viol_V36B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="36B"),])
  viol_V41A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="41A"),])
  viol_V42A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="42A"),])
  viol_V43A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="43A"),])
  viol_V44A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="44A"),])
  viol_V45A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="45A"),])
  viol_V46A <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="46A"),])
  viol_V41B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="41B"),])
  viol_V42B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="42B"),])
  viol_V43B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="43B"),])
  viol_V44B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="44B"),])
  viol_V45B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="45B"),])
  viol_V46B <- c(repl.raref[(substr(rownames(repl.raref),2,4)=="46B"),])
  # Draw boxplot-community
  pdf(paste(name_output,"_community.pdf",sep=""),width=12,height=12)
  color_names <- c("red3","orange","palegreen","cornflowerblue")
  par(fig=c(0,1,0.5,1))
  boxplot(c(viol_V11A,viol_V12A,viol_V13A,viol_V14A,viol_V16A),
          c(viol_V22A,viol_V23A,viol_V24A,viol_V25A,viol_V26A),
          c(viol_V31A,viol_V32A,viol_V33A,viol_V35A,viol_V36A),
          c(viol_V41A,viol_V42A,viol_V43A,viol_V44A,viol_V45A),
          col=grey.colors(4, start = 0.3, end = 1),ylab="MOTU richness (rarefied to 8000 reads)",
          names=c("V1A","V2A","V3A","V4A"),main="Fraction A (> 1 mm)",
          ylim=c(50,250))
  par(new=T,fig=c(0,1,0,0.5))
  boxplot(c(viol_V11B,viol_V12B,viol_V13B,viol_V14B,viol_V16B),
          c(viol_V22B,viol_V23B,viol_V24B,viol_V25B,viol_V26B),
          c(viol_V31B,viol_V32B,viol_V33B,viol_V35B,viol_V36B),
          c(viol_V41B,viol_V42B,viol_V43B,viol_V44B,viol_V45B),
          col=grey.colors(4, start = 0.3, end = 1),ylab="MOTU richness (rarefied to 8000 reads)",
          names=c("V1B","V2B","V3B","V4B"),main="Fraction B (63 µm - 1 mm)",
          ylim=c(50,250))
  dev.off()
  # Draw boxplot replicates
  pdf(paste(name_output,"_replicates.pdf",sep=""),width=15,height=12)
  par(fig=c(0,1,0.5,1))
  boxplot(viol_V11A,viol_V12A,viol_V13A,viol_V14A,viol_V15A,viol_V16A,
          viol_V21A,viol_V22A,viol_V23A,viol_V24A,viol_V25A,viol_V26A,
          viol_V31A,viol_V32A,viol_V33A,viol_V34A,viol_V35A,viol_V36A,
          viol_V41A,viol_V42A,viol_V43A,viol_V44A,viol_V45A,viol_V46A,
          col=c(rep(color_names[1],6),rep(color_names[3],6),rep(color_names[5],6),rep(color_names[7],6)),
          ylab="MOTU richness (rarefied to 8000 reads)",
          names=c(rep("V1A",6),rep("V2A",6),rep("V3A",6),rep("V4A",6)),main="Fraction A (> 1 mm)",
          ylim=c(50,250))
  par(new=T,fig=c(0,1,0,0.5))
  boxplot(viol_V11B,viol_V12B,viol_V13B,viol_V14B,viol_V15B,viol_V16B,
          viol_V21B,viol_V22B,viol_V23B,viol_V24B,viol_V25B,viol_V26B,
          viol_V31B,viol_V32B,viol_V33B,viol_V34B,viol_V35B,viol_V36B,
          viol_V41B,viol_V42B,viol_V43B,viol_V44B,viol_V45B,viol_V46B,
          col=c(rep(color_names[1],6),rep(color_names[3],6),rep(color_names[5],6),rep(color_names[7],6)),
          ylab="MOTU richness (rarefied to 8000 reads)",
          names=c(rep("V1B",6),rep("V2B",6),rep("V3B",6),rep("V4B",6)),main="Fraction B (63 µm - 1 mm)",
          ylim=c(50,250))
  dev.off()
  ### Test for differences
  replicates_list <- list(viol_V11A,viol_V12A,viol_V13A,viol_V14A,viol_V15A,viol_V16A,
                          viol_V21A,viol_V22A,viol_V23A,viol_V24A,viol_V25A,viol_V26A,
                          viol_V31A,viol_V32A,viol_V33A,viol_V34A,viol_V35A,viol_V36A,
                          viol_V41A,viol_V42A,viol_V43A,viol_V44A,viol_V45A,viol_V46A,
                          viol_V11B,viol_V12B,viol_V13B,viol_V14B,viol_V15B,viol_V16B,
                          viol_V21B,viol_V22B,viol_V23B,viol_V24B,viol_V25B,viol_V26B,
                          viol_V31B,viol_V32B,viol_V33B,viol_V34B,viol_V35B,viol_V36B,
                          viol_V41B,viol_V42B,viol_V43B,viol_V44B,viol_V45B,viol_V46B
  )
  medians <- unlist(lapply(replicates_list,median))
  medians_A <- medians[1:24]
  medians_B <- medians[25:48]
  factors <- as.factor(c(rep("V1A",6),rep("V2A",6),rep("V3A",6),rep("V4A",6),rep("V1B",6),rep("V2B",6),rep("V3B",6),rep("V4B",6)))
  factors_A <- factors[1:24]
  factors_B <- factors[25:48]
  dunn.test(medians_A,factors_A,method="bh",kw=T)
  dunn.test(medians_B,factors_B,method="bh",kw=T)
}

rarefaction_dunn_sr(db,"SWARM_VENT_all")

#### ANOVA_MDS ####
## Collapses species by name
  # Collapses species with the same name over an identity threshold
  # The input is a csv file with field names "best_identity", "species", "sequence", "total_reads".
  # Sample columns are located between start_samples and end_samples (numeric).
  # The output will have the most abundant sequence as the representative and the higher best_id as the best_id
  # Author: Owen S. Wangensteen, owenwangensteen@ub.edu
  # Metabarpark Project, metabarpark.blogspot.com
  # ChallenGen Project, challengenproject.blogspot.com

library("optparse")
db <- db_swarm_normalized
colnames(db)
names(db)[3] = "species"

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="csv file including 'id' and 'sequence' fields", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, 
              help="Output file name [default = input file ending in _collapsed.csv]", metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default= 0.70, 
              help="Threshold for collapsing. Default = 0.70", metavar="numeric"),
  make_option(c("-s", "--start_samples"), type="integer", default= 14, 
              help="Sample columns start. Default = 14", metavar="numeric"),
  make_option(c("-e", "--end_samples"), type="integer", default= 98, 
              help="Sample columns end. Default = 98", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile) ){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input csv file including 'best_identity', 'species', 'sequence' & 'total_reads' fields.", call.=FALSE)
}

# If no outfile name, then create a suitable one
if (is.null(opt$outfile)){
  opt$outfile <- paste(substr(opt$infile,1,nchar(opt$infile)-4),"_collapsed.csv",sep="")
}

# Read the csv file. Check that it has more than 1 column.
message("Reading ",opt$infile," database")
db <- na.exclude(read.table(opt$infile,sep=";",head=T))
if (ncol(db)==1) db <- na.exclude(read.table(opt$infile,sep=",",head=T))
if (ncol(db)==1) stop("Unable to read csv file. Use a file with ',' or ';' separators",call.=FALSE)
message("Database ",opt$infile," read including ",nrow(db), " sequences.")

# Check the number of species
unique_species <- length(unique(db$species))
db_species <- na.exclude(db[db$species!="",])
rows_with_species_name <- nrow(db_species)
message("Database ",opt$infile," includes ",rows_with_species_name, " sequences with species name,")
message("Belonging to ",unique_species," unique species.")
message("Now collapsing species with best_identity higher than ",opt$threshold,".")

# Select which sequences need to be collapsed
db_all_collapsed <- tapply(db_species$total_reads/db_species$total_reads,db_species$species,sum,na.rm=T)
db_to_collapse <- db[(db$species %in% names(db_all_collapsed[db_all_collapsed>1])) & (db$best_identity>opt$threshold),]
db_unmodified <-  db[!((db$species %in% names(db_all_collapsed[db_all_collapsed>1])) & (db$best_identity>opt$threshold)),]
species_to_collapse <- unique(as.character(db_to_collapse$species))

# Do the collapse calculations
species_collapsed <- NULL
for (name in species_to_collapse) {
  species_db <- db_to_collapse[db_to_collapse$species==name,]
  reads_collapsed <- colSums(species_db[,opt$start_samples:opt$end_samples])
  # Order dataframe by total_reads and best_id
  species_db <- species_db[with(species_db,order(total_reads,best_identity,decreasing=T)),]
  # Get the first row
  vector <- species_db[1,]
  # Change sample columns with the sum
  vector[1,opt$start_samples:opt$end_samples] <- reads_collapsed
  # Change total_reads
  vector$total_reads[1] <- sum(reads_collapsed)
  # Change best_id
  vector$best_identity[1] <- max(species_db$best_identity)
  species_collapsed <- rbind(species_collapsed,vector)
}
message(nrow(db_to_collapse)," sequences collapsed to ",length(species_to_collapse)," species.")

# Recover the database
db_final <- rbind(species_collapsed,db_unmodified)
db_final <- db_final[order(db_final$superkingdom,db_final$kingdom,db_final$phylum,db_final$class,db_final$order,
                           db_final$family,db_final$genus,db_final$species,db_final$scientific_name,db_final$group1,db_final$group2, db_final$group3,db_final$total_reads),]
write.table(db_final,file=opt$outfile,append = F,row.names=F,sep=",")
message("Written ",opt$outfile," with ",nrow(db_final)," MOTUs.")

# END #

### ANOVA_ALL 
db <- read.table("~/Desktop/metabarcodinganalysis/VENT_SWARM_collapsed.csv",header = T, sep = ",", stringsAsFactors = FALSE)
colnames(db)

## Function for pairwise-adonis
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m, permut=9999)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~ factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] ,
                method =sim.method,
                permutations=permut)
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

# ANOVA_analysis

sample_columns <- c(17:64)
abundances <- db[,sample_columns]
#PERANOVA(ADONIS) #four-root and Brya curtis
total_abundances <- rowSums(abundances)
gradient <- substr(colnames(abundances),1,2)
fraction <- substr(colnames(abundances),4,4)
gradient <- as.factor(gradient)
fraction <- as.factor(fraction)

abundances_fraction <- abundances[,fraction=="A"]
abundances_corrected <- abundances_fraction[rowSums(abundances_fraction)>0,]
results <- adonis(t(((sqrt(sqrt(renormalize(abundances_corrected))))))~gradient[fraction=="A"],method="bray",permutations = 9999)
results
# Calculate pairwise-adonis and write the output table to a file
PW.Adonis_results=pairwise.adonis(x = t(((((renormalize(abundances_corrected)))))),factors = gradient[fraction=="A"],sim.method="bray",p.adjust.m = "BH",permut = 9999)
write.table(PW.Adonis_results,"pairwise-adonis_allA.csv",sep=",",row.names=F) 

abundances_fraction <- abundances[,fraction=="B"]
abundances_corrected <- abundances_fraction[rowSums(abundances_fraction)>0,]
results <- adonis(t(((sqrt(sqrt(renormalize(abundances_corrected))))))~gradient[fraction=="B"],method="bray",permutations = 9999)
results
# Calculate pairwise-adonis and write the output table to a file
PW.Adonis_results=pairwise.adonis(x = t(((((renormalize(abundances_corrected)))))),factors = gradient[fraction=="B"],sim.method="bray",p.adjust.m = "BH",permut = 9999)
write.table(PW.Adonis_results,"pairwise-adonis_allB.csv",sep=",",row.names=F) 

### MDS 
library(ade4)
library(vegan)
library(MASS)

db <- read.table("~/Desktop/metabarcodinganalysis/VENT_SWARM_collapsed.csv",header = T, sep = ",", stringsAsFactors = FALSE)

colnames(db)
sample_columns <- 17:64
db$group1

#1-Algae#
group1="Algae"

# FRACTION A
# MDS
abundances <- db[tolower(db$group1)==group,sample_columns]
row.names(abundances) <- db$id[tolower(db$group1)==group]

# FRACTION A
total_abundances <- rowSums(abundances)
gradient <- substr(colnames(abundances),1,2)
fraction <- substr(colnames(abundances),4,4)
gradient <- as.factor(gradient)
fraction <- as.factor(fraction)
abundances_fractionA <- abundances[,fraction=="A"]
abundances_correctedA <- abundances_fractionA[rowSums(abundances_fractionA)>0,]

distancesA <- vegdist(t(((sqrt(sqrt(renormalize(abundances_correctedA)))))),method = "bray")
mdsA <- isoMDS(distancesA,maxit = 300)

# PLOT MDS #
#mds$points[,1] <- -mds$points[,1]
draw_mdsA <- function(code,codesgroup,col,pch,cex,label){
  points(mdsA$points[code,], pch=pch,col=col,cex=cex)
  text(mdsA$points[code,],labels = label,cex=.8)
}

codes <- as.factor(paste0(substr(colnames(abundances_correctedA),1,2),"_",substr(colnames(abundances_correctedA),4,5)))
V1A <- codes=="V1_A"
V2A <- codes=="V2_A"
V3A <- codes=="V3_A"
V4A <- codes=="V4_A"
pdf("MDS_allA.pdf",width=12,height=12)
par(fig=c(0,1,0,1),mar=c(5,5,5,5))
color_names <- c( "black", "red3","orange","palegreen","cornflowerblue")
plot(mdsA$points, type="n",xaxt="n",yaxt="n",xlab="MDS 1",ylab="MDS 2", cex.lab=2,
     xlim=c(min(mdsA$points[,1])-.2,max(mdsA$points[,1])+.2),
     ylim=c(min(mdsA$points[,2])-.2,max(mdsA$points[,2])+.2),
     main=paste0("MDS A. ",nrow(abundances_correctedA)," MOTUs"))
draw_mdsA(V1A,"V1_A",color_names[1],15,2,"")
draw_mdsA(V2A,"V2_A",color_names[1],8,2,"")
draw_mdsA(V3A,"V3_A",color_names[1],3,2,"")
draw_mdsA(V4A,"V4_A",color_names[1],0,2,"")

legend("bottomleft",legend=c("1-Vent","2-Transition25","3-Transition75","4-Control"),pch=c(15,8,3,0),cex=2,bty="o")
legend("top",legend=c(paste("stress = ",round(mdsA$stress,2)," %",sep="")),bty="n",cex=2)
plot(efA, p.max = 0.05, col = "black" )
dev.off()

# FRACTION B
abundances_fractionB <- abundances[,fraction=="B"]
abundances_correctedB <- abundances_fractionB[rowSums(abundances_fractionB)>0,]

distancesB <- vegdist(t(((sqrt(sqrt(renormalize(abundances_correctedB)))))),method = "bray")
mdsB <- isoMDS(distancesB,maxit = 300)

# PLOT MDS #
#mds$points[,1] <- -mds$points[,1]
draw_mdsB <- function(code,codesgroup,col,pch,cex,label){
  points(mdsB$points[code,], pch=pch,col=col,cex=cex)
  text(mdsB$points[code,],labels = label,cex=.8)
}

codes <- as.factor(paste0(substr(colnames(abundances_correctedB),1,2),"_",substr(colnames(abundances_correctedB),4,5)))
V1B <- codes=="V1_B"
V2B <- codes=="V2_B"
V3B <- codes=="V3_B"
V4B <- codes=="V4_B"
pdf("MDS_allB.pdf",width=12,height=12)
par(fig=c(0,1,0,1),mar=c(5,5,5,5))
color_names <- c("black", "red3","orange","palegreen","cornflowerblue")
plot(mdsB$points, type="n",xaxt="n",yaxt="n",xlab="MDS 1",ylab="MDS 2", cex.lab=2,
     xlim=c(min(mdsB$points[,1])-.2,max(mdsB$points[,1])+.2),
     ylim=c(min(mdsB$points[,2])-.2,max(mdsB$points[,2])+.2),
     main=paste0("MDS B. ",nrow(abundances_correctedB)," MOTUs"))
draw_mdsB(V1B,"V1_B",color_names[1],15,2,"")
draw_mdsB(V2B,"V2_B",color_names[1],8,2,"")
draw_mdsB(V3B,"V3_B",color_names[1],3,2,"")
draw_mdsB(V4B,"V4_B",color_names[1],0,2,"")

legend("bottomleft",legend=c("1-Vent","2-Transition25","3-Transition75","4-Control"),pch=c(15,8,3,0),cex=2,bty="o")
legend("top",legend=c(paste("stress = ",round(mdsB$stress,2)," %",sep="")),bty="n",cex=2)
plot(efB, p.max = 0.05, col = "black" )
dev.off()

#### SIMPER ####
db <- read.table("~/Desktop/metabarcodinganalysis/VENT_SWARM_collapsed.csv",header = T, sep = ",", stringsAsFactors = FALSE)
colnames(db)
sample_columns <- c(17:64)
abundances <- db[,sample_columns]
row.names(abundances) <- db$species
gradient <- substr(colnames(abundances),1,2)
fraction <- substr(colnames(abundances),4,4)
gradient <- as.factor(gradient)
fraction <- as.factor(fraction)

# Fraction A
abundances_fraction <- abundances[,fraction=="A"]
abundances_corrected <- abundances_fraction[rowSums(abundances_fraction)>0,]
com_matrix <- t(sqrt(sqrt(renormalize(abundances_corrected))))
codes <- gradient
simper_AllSamples <- simper(com_matrix,gradient[fraction=="A"])
# V1_V2
# with 30 species per comparision
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V1_V2"]]$species[as.numeric(simper_AllSamples[["V1_V2"]]$ord[1:30])]
for (i in 1:30) table_sum$rank[i] <- as.character(db$rank[db$id==table_sum$id[i]])
for (i in 1:30) table_sum$scientific_name[i] <- as.character(db$species[db$id==table_sum$id[i]])
for (i in 1:30) table_sum$best_id[i] <- (db$best_identity[db$id==table_sum$id[i]])
for (i in 1:30) table_sum$superkingdom[i] <- as.character(db$superkingdom[db$id==table_sum$id[i]])
for (i in 1:30) table_sum$kingdom[i] <- as.character(db$kingdom[db$id==table_sum$id[i]])
for (i in 1:30) table_sum$group[i] <- as.character(db$group[db$id==table_sum$id[i]])
table_sum$ava <- round(100*simper_AllSamples[["V1_V2"]]$ava[simper_AllSamples[["V1_V2"]]$species[as.numeric(simper_AllSamples[["V1_V2"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V1_V2"]]$avb[simper_AllSamples[["V1_V2"]]$species[as.numeric(simper_AllSamples[["V1_V2"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V1_V2"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V1_V2A.csv",row.names=F,sep=",",quote = F)
#V1_V3
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V1_V3"]]$species[as.numeric(simper_AllSamples[["V1_V3"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V1_V3"]]$ava[simper_AllSamples[["V1_V3"]]$species[as.numeric(simper_AllSamples[["V1_V3"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V1_V3"]]$avb[simper_AllSamples[["V1_V3"]]$species[as.numeric(simper_AllSamples[["V1_V3"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V1_V3"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V1_V3A.csv",row.names=F,sep=",",quote = F)
#V1_V4
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V1_V4"]]$species[as.numeric(simper_AllSamples[["V1_V4"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V1_V4"]]$ava[simper_AllSamples[["V1_V4"]]$species[as.numeric(simper_AllSamples[["V1_V4"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V1_V4"]]$avb[simper_AllSamples[["V1_V4"]]$species[as.numeric(simper_AllSamples[["V1_V4"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V1_V4"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V1_V4A.csv",row.names=F,sep=",",quote = F)
#V2_V3
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V2_V3"]]$species[as.numeric(simper_AllSamples[["V2_V3"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V2_V3"]]$ava[simper_AllSamples[["V2_V3"]]$species[as.numeric(simper_AllSamples[["V2_V3"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V2_V3"]]$avb[simper_AllSamples[["V2_V3"]]$species[as.numeric(simper_AllSamples[["V2_V3"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V2_V3"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V2_V3A.csv",row.names=F,sep=",",quote = F)
#V2_V4
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V2_V4"]]$species[as.numeric(simper_AllSamples[["V2_V4"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V2_V4"]]$ava[simper_AllSamples[["V2_V4"]]$species[as.numeric(simper_AllSamples[["V2_V4"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V2_V4"]]$avb[simper_AllSamples[["V2_V4"]]$species[as.numeric(simper_AllSamples[["V2_V4"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V2_V4"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V2_V4A.csv",row.names=F,sep=",",quote = F)
#V3_V4
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V3_V4"]]$species[as.numeric(simper_AllSamples[["V3_V4"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V3_V4"]]$ava[simper_AllSamples[["V3_V4"]]$species[as.numeric(simper_AllSamples[["V3_V4"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V3_V4"]]$avb[simper_AllSamples[["V3_V4"]]$species[as.numeric(simper_AllSamples[["V3_V4"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V3_V4"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V3_V4A.csv",row.names=F,sep=",",quote = F)

# Fraction B
abundances_fraction <- abundances[,fraction=="B"]
abundances_corrected <- abundances_fraction[rowSums(abundances_fraction)>0,]
com_matrix <- t(sqrt(sqrt(renormalize(abundances_corrected))))
codes <- gradient
simper_AllSamples <- simper(com_matrix,gradient[fraction=="B"])
# V1_V2
# with 30 species per comparision
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V1_V2"]]$species[as.numeric(simper_AllSamples[["V1_V2"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V1_V2"]]$ava[simper_AllSamples[["V1_V2"]]$species[as.numeric(simper_AllSamples[["V1_V2"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V1_V2"]]$avb[simper_AllSamples[["V1_V2"]]$species[as.numeric(simper_AllSamples[["V1_V2"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V1_V2"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V1_V2B.csv",row.names=F,sep=",",quote = F)
#V1_V3
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V1_V3"]]$species[as.numeric(simper_AllSamples[["V1_V3"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V1_V3"]]$ava[simper_AllSamples[["V1_V3"]]$species[as.numeric(simper_AllSamples[["V1_V3"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V1_V3"]]$avb[simper_AllSamples[["V1_V3"]]$species[as.numeric(simper_AllSamples[["V1_V3"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V1_V3"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V1_V3B.csv",row.names=F,sep=",",quote = F)
#V1_V4
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V1_V4"]]$species[as.numeric(simper_AllSamples[["V1_V4"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V1_V4"]]$ava[simper_AllSamples[["V1_V4"]]$species[as.numeric(simper_AllSamples[["V1_V4"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V1_V4"]]$avb[simper_AllSamples[["V1_V4"]]$species[as.numeric(simper_AllSamples[["V1_V4"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V1_V4"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V1_V4Bcsv",row.names=F,sep=",",quote = F)
#V2_V3
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V2_V3"]]$species[as.numeric(simper_AllSamples[["V2_V3"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V2_V3"]]$ava[simper_AllSamples[["V2_V3"]]$species[as.numeric(simper_AllSamples[["V2_V3"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V2_V3"]]$avb[simper_AllSamples[["V2_V3"]]$species[as.numeric(simper_AllSamples[["V2_V3"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V2_V3"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V2_V3B.csv",row.names=F,sep=",",quote = F)
#V2_V4
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V2_V4"]]$species[as.numeric(simper_AllSamples[["V2_V4"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V2_V4"]]$ava[simper_AllSamples[["V2_V4"]]$species[as.numeric(simper_AllSamples[["V2_V4"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V2_V4"]]$avb[simper_AllSamples[["V2_V4"]]$species[as.numeric(simper_AllSamples[["V2_V4"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V2_V4"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V2_V4B.csv",row.names=F,sep=",",quote = F)
#V3_V4
table_sum <- data.frame(id = rep("",30),rank= rep("",30),scientific_name = rep("",30),best_id = rep(0,30),superkingdom = rep("",30),kingdom = rep("",30),group = rep("",30),ava= rep(0,30),avb= rep(0,30),cumsum= rep(0,30),stringsAsFactors = F)
table_sum$id <- simper_AllSamples[["V3_V4"]]$species[as.numeric(simper_AllSamples[["V3_V4"]]$ord[1:30])]
table_sum$ava <- round(100*simper_AllSamples[["V3_V4"]]$ava[simper_AllSamples[["V3_V4"]]$species[as.numeric(simper_AllSamples[["V3_V4"]]$ord[1:30])]],3)
table_sum$avb <- round(100*simper_AllSamples[["V3_V4"]]$avb[simper_AllSamples[["V3_V4"]]$species[as.numeric(simper_AllSamples[["V3_V4"]]$ord[1:30])]],3)
table_sum$cumsum <- round(100*simper_AllSamples[["V3_V4"]]$cusum[1:30],3)
write.table(table_sum,"Simper_V3_V4B.csv",row.names=F,sep=",",quote = F)

#### Barplot ####
db <- read.table("~/Desktop/metabarcodinganalysis/VENT_SWARM_collapsed.csv",header = T, sep = ",", stringsAsFactors = FALSE)
colnames(db)
#Change the order of the columns (Replicates) 
library(dplyr)
db2 <- db %>% select(group3,V11A,V12A,V13A,V14A,V15A,V16A,V11B,V12B,V13B,V14B,V15B,V16B,V21A,
                     V22A,V23A,V24A,V25A,V26A,V21B,V22B,V23B,V24B,V25B,V26B,V31A,V32A,V33A,V34A,
                     V35A,V36A,V31B,V32B,V33B,V34B,V35B,V36B,V41A,V42A,V43A,V44A,V45A,V46A,V41B,
                     V42B,V43B,V44B,V45B,V46B)
colnames(db2)
#Replicates sums
sample_columns <- c(2:49)
group3 <-db2[,1] 
db2_sum <- data.frame(matrix(0,nrow=8,ncol=length(seq(2,49,6)))) 
for (i in 1:8) db2_sum[,i] <- rowSums(db2[,((2+6*(i-1)):(7+6*(i-1)))])
namescol <- NULL
for (i in 1:8) namescol <- c(namescol,substr(colnames(db2[(2+6*(i-1))]),2,4)) 
colnames(db2_sum) <- namescol
db2_sum <- cbind(group3,db2_sum)

#### BUBBLE_plot ####
#Select the sp of interest from VENT_SWARM_collapsed.csv
#Obtain the values of MDS dimensions 1 and 2 to created the bubble plot.   
library(plotly)
df -> read.table("~/Desktop/metabarcodinganalysis/VENT_Eurythoe.csv",header = T, sep = ",", stringsAsFactors = FALSE)
bubbleplot <- plot_ly(df, x = ~mds1, y = ~mds2,
                      text = ~Fraction, size = ~abundances,
                      color = ~Fraction, color = I("brown3","brown1"),
                      sizes = c(10, 50),
                      marker =
                        list(opacity = 0.5,
                             sizemode = "diameter"))

bubbleplot <- bubbleplot%>%layout (title = "Eurythoe", xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FAlSE))
bubbleplot

#### HEATMAP ####
#The input document is a csv with means reads per sample
#it has to be a square matrix
library(gplots)
library(ggthemes) 
library(readr)

y <- read.table("~/Desktop/metabarcodinganalysis/matrix_heatmap.csv",header = T, sep = ",", stringsAsFactors = FALSE)

colfunc<- colorRampPalette(c("white","black"))
colfunc(4)
my_palette<-c("#FFFFFF", "#AAAAAA", "#555555", "#000000")
my_palette
heatmap.2(y,Rowv=TRUE, Colv=TRUE, legend=TRUE, dendrogram=("both"), col=my_palette, density.info="none", trace="none")
?heatmap.2
