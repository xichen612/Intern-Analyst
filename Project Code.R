# Name
# ====
# Xi Chen
# Statistics, MS
# Univeristy of California, Davis

# Project, Title
# ==============
# Intern. RA, Dept. of Ecosystem Science and Management, PSU

# Instructor
# ==========
# Prof.Henry Lin

# Time
# ====
# 2015.6.15-2015.8.15.

# Content
# =======
# Data Anlysis on three Soil Databases (ISRIC-WISE, HWSD, NRCS)
#   * Clean Database, Combing Tables, Droping Outliers.
#   * Histogram and Density Plot of Global Distribution of Soil Propertities in each database.
#   * Spearman Correlation of all properties in each database.
#   * PTFs: Bulk.Density and CEC.Soil(Clay) Prediction using other properties.
#   * Clustering by all Chemical and Physical Properties.
#   * Soil Set, Soil Group, Soil subset Distribution.
#   * Soil Set(subset) Classification using their propertites

# Work Hard and Be Creative.

# Functions.

###Clean data###
# grap the common name of each column
subnames = function(x){
  # Args:
  # * x: colnames of dataset
  
  for(i in 1:length(x)){
    x[i] = gsub("^[[:upper:]]_", "", colnames(sep1)[i])
  }
  return(x)
}

#extract hdat-related set from the SMU subset 
smubind = function(x){
  # Args:
  # * x: data
  
  cols = dim(x)[1]
  m = NULL
  for(i in 1:cols){
    y = smu[which(smu[,2] == x[i,2]),]
    m = rbind(m, y)
  }
  return(m)
}


###global distribution### 
# histogram with normal curve
# example:
# normHist(temp, 2, 4.5,4.3,4.1, .05,un6, 7, 12,15,"black")
library(ggplot2)
normHist = function(dat, trans, a, b1, b2, b3, width, un, ss, s1, s2, colors){
  # Args:
  # * dat: clean data (eg: Horizon and Sand)
  # * trans: what kind of transformation to the dat
  # * a: Horizontal label.
  # * b1: Vertical coordinate of mean.
  # * b2: Vertical coordinate of sd.
  # * width: width of histogram
  # * un: unit of property.
  # * ss: size of note.(mean, sd, N)
  # * s1: size of axis
  # * s2: size of axis title
  # * colors: color of note.
  
  m = round(mean(dat[,2]), digits = 3)
  s = round(sd(dat[,2]), digits = 3)
  n = dim(dat)[1]
  unit = un
  label1 = paste("mean =", m)
  label2 = paste("sd =", s)
  label3 = paste("N =", n)
  input = paste(colnames(dat)[2], unit)
  colnames(dat)[2] = "Property"
  
  if(trans == "log"){
    dat[,2] = log10(dat[,2])
    input = paste("Log", input)
  }else{
    if(trans == "sqrt"){
      dat[,2] = sqrt(dat[,2])
      input = paste("Sqrt", input)
    }
  }
  me1 = round(mean(dat[,2]), digits = 3)
  sd1 = round(sd(dat[,2]), digits = 3)
  
  gg = ggplot(dat, aes(x = Property)) + 
    geom_histogram(aes(y = ..density..,  fill=..count..), binwidth=width) + 
    scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C") +
    stat_function(fun = dnorm, colour = "red",  
                  args = list(mean = me1, sd = sd1)) +
    ylab("Relative Frequency") +
    xlab(input) +
    annotate("text", x = c(rep(a,3)), y = c(b1, b2, b3), label = c(label1, label2, label3), 
             size = ss, colour = colors) +
    theme(axis.text=element_text(size=s1),axis.title=element_text(size=s2,face="bold"))
  
  return(gg)
}


# boxplot with Tukey pairwise comparison
library(agricolae)
# for A/B/C Horizon
TukBox = function(dat, unit, v, s, t1, t2, trans){
    # Args:
    # * dat: data only includes property, Order and Horizon by first, second and third column.(no NA)
    # * unit: the unit of property
    # * v: vertical coordinate of Tukey's comparison label.
    # * s: size of the label.
    # * t1, t2: size of axis text and axis title.
    
    input = paste(colnames(dat)[2], unit)
    colnames(dat)[2] = "Property"
    model = aov(Property~Horizon, data = dat)
    tuk = HSD.test(model, "Horizon", group = TRUE, alpha = .05)
    labels = tuk$group[,3]
    len = length(labels)
    tables = table(dat[,1])
    ho = as.character(dat[,1])
    ho[which(ho == "A")] = paste0("A (n = ",tables[1], ")")
    ho[which(ho == "B")] = paste0("B (n = ",tables[2], ")")
    ho[which(ho == "C")] = paste0("C (n = ",tables[3], ")")
    ho = as.factor(ho)
    dat[,1] = ho
    gg = ggplot(dat, aes(x=Horizon, y=Property, fill= Horizon)) + geom_boxplot() +
      ylab(input)+
      xlab("")+
      annotate("text", x = c(1:len), y = rep(v, len), label = labels, 
               size = s, colour = "black") +
      theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
    
    if(trans == "log"){
      gg = gg + scale_y_continuous(trans=log10_trans())
    }
    
    return(gg)
  }  
#eg: TukBox(temp, un5, 105, 8, 15, 25)
# for Top/Sub Horizon
TukBox2 = function(dat, unit, v, s, t1, t2, trans){
  # Args:
  # * dat: data only includes property, Order and Horizon by first, second and third column.(no NA)
  # * unit: the unit of property
  # * v: vertical coordinate of Tukey's comparison label.
  # * s: size of the label.
  # * t1, t2: size of axis text and axis title.
  
  input = paste(colnames(dat)[2], unit)
  colnames(dat)[2] = "Property"
  model = aov(Property~Horizon, data = dat)
  tuk = HSD.test(model, "Horizon", group = TRUE, alpha = .05)
  labels = tuk$group[,3]
  len = length(labels)
  tables = table(dat[,1])
  ho = as.character(dat[,1])
  ho[which(ho == "Top")] = paste0("Top (n = ",tables[1], ")")
  ho[which(ho == "Sub")] = paste0("Sub (n = ",tables[2], ")")
  ho = as.factor(ho)
  dat[,1] = ho
  dat$Horizon <- factor(dat$Horizon, levels = rev(levels(dat$Horizon)))
  gg = ggplot(dat, aes(x=Horizon, y=Property, fill= Horizon)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
    annotate("text", x = c(1:len), y = rep(v, len), label = labels, 
             size = s, colour = "black") +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  
  if(trans == "log"){
    gg = gg + scale_y_continuous(trans=log10_trans())
  }else{
    if(trans == "sqrt"){
      gg = gg + scale_y_continuous(trans=sqrt_trans())
    }
  }
  
  return(gg)
} 

TukeyBoxplot = function(dat, x, c1, c2, c3, c4, a1, a2, a3, a4, s1, s2, unit){
  # Args:
  # * x: Vertical position.
  # * c1, c2, c3, c4: colors of mean, max, min, pvalue.
  # * a1, a2, a3, a4: size of mean, max, min, pvalue
  # * s1, s2: size of axis and axis title.
  # * unit: unit of property.
  
  input = paste(colnames(dat)[2], unit)
  colnames(dat)[2] = "Property" 
  colnames(dat)[1] = "Horizon"
  model <- aov(Property ~ Horizon, data = dat)
  tt = TukeyHSD(model)$Horizon[,4]
  nn = names(tt)
  tt = round(tt, digits = 3)
  p = NULL
  g = ": Pvlaue = "
  for(i in 1:length(tt)){
    p[i] = paste0(nn[i], g)
    p[i] = paste0(p[i], tt[i])
  }
  
  tables = table(dat[,1])
  ho = as.character(dat[,1])
  ho[which(ho == "A")] = paste0("A (n = ",tables[1], ")")
  ho[which(ho == "B")] = paste0("B (n = ",tables[2], ")")
  ho[which(ho == "C")] = paste0("C (n = ",tables[3], ")")
  ho = as.factor(ho)
  dat[,1] = ho  
  
  gg = ggplot(dat, aes(x=Horizon, y=Property, fill=Horizon)) + geom_boxplot()+ 
    ylab(input) +
    xlab("") +
    stat_summary(fun.y= mean, colour="darkred", geom="point", shape=18, size=3, show_guide = FALSE) +
    stat_summary(fun.y= mean, colour=c1, geom="text", size = a1, show_guide = FALSE, 
                 vjust=-0.7, aes( label=round(..y.., digits=1))) +
    stat_summary(fun.y=max, colour=c2, geom="text", size = a2, show_guide = FALSE, 
                 vjust=-0.2, aes( label=round(..y.., digits=1))) +
    stat_summary(fun.y=min, colour=c3, geom="text", size = a3, show_guide = FALSE, 
                 vjust=-0.7, aes( label=round(..y.., digits=1))) +
    annotate("text", x = c(1:3), y = rep(x,3), label = p, colour = c4, size = a4) +
    theme(axis.text=element_text(size=s1),axis.title=element_text(size=s2,face="bold"))
  
  return(gg)
}



###the correlation plot###
# put histogram on the diagnal
panel.hist1 <- function(x, ...)
{ usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
# put correlations on the upper panels,
panel.cor1 <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs", method = "spearman")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}
# regression on the lower panels.
panel.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex) 
  ok <- is.finite(x) & is.finite(y) 
  if (any(ok)) 
    abline(stats::lm(y[ok] ~ x[ok]), col = col.regres, ...) 
}
# eg
# pairs(c1, lower.panel=panel.smooth, upper.panel=panel.cor1, diag.panel=panel.hist1)



###SoilOrders comparison###
compBoxplot(temp,1,105,3.5, 5, 12,25,"black", un5)
compBoxplot = function(d,a,b1,b2,s,t1,t2,colors,sen){
  # Args:
  # * a: horizon axis of note(p value)
  # * b1, b2: vertiacal axis of note ( seq(b1, b1+(len-1)*b2, by = b2))
  # * s: size of note.
  # * colors: color of note.
  # * t1, t2: size of axis and axis title.
  
  name = colnames(d)[2]
  input = paste(name, sen)
  colnames(d)[2] = "Property"
  ttt = aov(Property~Order, data = d)
  if(all(TukeyHSD(ttt)$Order[,4] < 0.05)){
    gg = ggplot(d, aes(x=Order, y=Property, fill=Horizon)) + geom_boxplot() +
      ylab(input)+
      xlab("")+
      theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  }else{
  re = TukeyHSD(ttt)$Order[TukeyHSD(ttt)$Order[,4] > 0.05, 4]
  ttname = row.names(TukeyHSD(ttt)$Order)[TukeyHSD(ttt)$Order[,4] > 0.05]
  re = round(re, digits = 3)
  len = length(ttname)
  p = NULL
  for(i in 1:len){
    p[i] = paste(ttname[i], ":",re[i])
  }
  gg = ggplot(d, aes(x=Order, y=Property, fill=Horizon)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
    annotate("text", x = c(rep(a,len)), y = c(seq(b1, b1+(len-1)*b2, by = b2)), label = p, 
             size = s, colour = colors) +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  }
  return(gg)
}

# this is for Top/Sub Horizon database
ordercomp = function(dat, unit, b, c, s, t1, t2){
  # Args:
  # * dat: data only includes property, horizon and order.
  # * unit: the unit of property
  # * b: vertical coordinate of Tukey's comparison label.of (Top Horizon)
  # * c: vertical coordinate of Tukey's comparison label.of (Sub Horizon)
  # * s: size of the label.
  # * t1, t2: size of axis text and axis title.
  input = paste0(colnames(dat)[1], unit)
  colnames(dat)[1] = "Property"
  Topdat = dat[which(dat[,3] == "Top"), ]
  Subdat = dat[which(dat[,3] == "Sub"), ]
  model1 = aov(Property~Order, data = Topdat)
  tuk1 = HSD.test(model1, "Order", group=TRUE, alpha = 0.05)
  model2 = aov(Property~Order, data = Subdat)
  tuk2 = HSD.test(model2, "Order", group=TRUE, alpha = 0.05)
  labels1 = tuk1$group[,3]
  labels2 = tuk2$group[,3]
  len = length(levels(dat$Order))
  gg = ggplot(dat, aes(x=Order, y=Property, fill=Horizon)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
    annotate("text", x = c(1:len), y = rep(b, len), label = labels1, 
             size = s, colour = "black") +
    annotate("text", x = c(1:len), y = rep(c, len), label = labels2, 
             size = s, colour = "black") +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  
  return(gg)
}
# this is for A/B/C Horizon database
ordercomp2 = function(dat, unit, b, c, s, t1, t2, v1, v2, v3, w){
  # Args:
  # * dat: data only includes property, Order and Horizon by first, second and third column.(no NA)
  # * unit: the unit of property
  # * b: vertical coordinate of Tukey's comparison label.of (Top Horizon)
  # * c: vertical coordinate of Tukey's comparison label.of (Sub Horizon)
  # * s: size of the label.
  # * t1, t2: size of axis text and axis title.
  # * v1, v2, v3, w: adjustment of position of label(group letter) of Tukey's comparison.
  
  input = paste0(colnames(dat)[1], unit)
  colnames(dat)[1] = "Property"
  Adat = dat[which(dat[,3] == "A"), ]
  Bdat = dat[which(dat[,3] == "B"), ]
  Cdat = dat[which(dat[,3] == "C"), ]
  model1 = aov(Property~Order, data = Adat)
  tuk1 = HSD.test(model1, "Order", group=TRUE, alpha = 0.05)
  model2 = aov(Property~Order, data = Bdat)
  tuk2 = HSD.test(model2, "Order", group=TRUE, alpha = 0.05)
  model3 = aov(Property~Order, data = Cdat)
  tuk3 = HSD.test(model3, "Order", group=TRUE, alpha = 0.05)
  labels1 = tuk1$group[,3]
  labels2 = tuk2$group[,3]
  labels3 = tuk3$group[,3]
  len = length(levels(dat$Order))
  gg = ggplot(dat, aes(x=Horizon, y=Property, fill= Order)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
    #    annotate("text", x = c(1:len), y = rep(b, len), label = labels1, 
    #             size = s, colour = "black") +
    #    annotate("text", x = c(1:len), y = rep(c, len), label = labels2, 
    #             size = s, colour = "black") +
    annotate("text", x = c(seq(v1, v1+(len-1)*w, by = w)), y = rep(c, len), label = labels1, 
             size = s, colour = "black") +
    annotate("text", x = c(seq(v2, v2+(len-1)*w, by = w)), y = rep(b, len), label = labels2, 
             size = s, colour = "black") +
    annotate("text", x = c(seq(v3, v3+(len-1)*w, by = w)), y = rep(c, len), label = labels3, 
             size = s, colour = "black") +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  
  return(gg)
}
# eg: NRCS database sand:
# temp = na.omit(dat2[,c(2,11,1)]); head(temp)
# temp = temp[which(temp[,1]<100),]
# ordercomp2(temp,un5, 105,-5, 5, 12,25, 0.655, 1.655, 2.655, 0.063)
ordercomp3 = function(dat, unit, b, c, d, s, t1, t2){
  # Args:
  # * dat: data only includes property, Order and Horizon by first, second and third column.(no NA)
  # * unit: the unit of property
  # * b: vertical coordinate of Tukey's comparison label.of (A Horizon)
  # * c: vertical coordinate of Tukey's comparison label.of (B Horizon)
  # * d: vertical coordinate of Tukey's comparison label.of (C Horizon)
  # * s: size of the label.
  # * t1, t2: size of axis text and axis title.
  
  input = paste0(colnames(dat)[1], unit)
  colnames(dat)[1] = "Property"
  Adat = dat[which(dat[,3] == "A"), ]
  Bdat = dat[which(dat[,3] == "B"), ]
  Cdat = dat[which(dat[,3] == "C"), ]
  model1 = aov(Property~Order, data = Adat)
  tuk1 = HSD.test(model1, "Order", group=TRUE, alpha = 0.05)
  model2 = aov(Property~Order, data = Bdat)
  tuk2 = HSD.test(model2, "Order", group=TRUE, alpha = 0.05)
  model3 = aov(Property~Order, data = Cdat)
  tuk3 = HSD.test(model3, "Order", group=TRUE, alpha = 0.05)
  labels1 = tuk1$group[,3]
  labels2 = tuk2$group[,3]
  labels3 = tuk3$group[,3]
  len1 = length(labels1)
  len2 = length(labels2)
  len3 = length(labels3)
  gg = ggplot(dat, aes(x=Order, y=Property, fill= Horizon)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
        annotate("text", x = c(1:len1), y = rep(b, len1), label = labels1, 
                 size = s, colour = "black") +
        annotate("text", x = c(1:len2), y = rep(c, len2), label = labels2, 
                 size = s, colour = "black") +
        annotate("text", x = c(1:len3), y = rep(d, len3), label = labels3, 
                 size = s, colour = "black") +
#    annotate("text", x = c(seq(v1, v1+(len-1)*w, by = w)), y = rep(c, len), label = labels1, 
#             size = s, colour = "black") +
#    annotate("text", x = c(seq(v2, v2+(len-1)*w, by = w)), y = rep(b, len), label = labels2, 
#             size = s, colour = "black") +
#    annotate("text", x = c(seq(v3, v3+(len-1)*w, by = w)), y = rep(c, len), label = labels3, 
#             size = s, colour = "black") +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  
  return(gg)
}

###Soil set comparison###
# For Horizon A/B/C
setcomp = function(dat, unit, b, c, d, s, t1, t2, trans){
  # Args:
  # * dat: data only includes property, Order and Horizon by first, second and third column.(no NA)
  # * unit: the unit of property
  # * b: vertical coordinate of Tukey's comparison label.of (A Horizon)
  # * c: vertical coordinate of Tukey's comparison label.of (B Horizon)
  # * d: vertical coordinate of Tukey's comparison label.of (C Horizon)
  # * s: size of the label.
  # * t1, t2: size of axis text and axis title.
  # * trans: what transformation to the data, "log" stands for log10 transformation.
  
  input = paste0(colnames(dat)[1], unit)
  colnames(dat)[1] = "Property"
  Adat = dat[which(dat[,3] == "A"), ]
  Bdat = dat[which(dat[,3] == "B"), ]
  Cdat = dat[which(dat[,3] == "C"), ]
  model1 = aov(Property~Order, data = Adat)
  tuk1 = HSD.test(model1, "Order", group=TRUE, alpha = 0.05)
  model2 = aov(Property~Order, data = Bdat)
  tuk2 = HSD.test(model2, "Order", group=TRUE, alpha = 0.05)
  model3 = aov(Property~Order, data = Cdat)
  tuk3 = HSD.test(model3, "Order", group=TRUE, alpha = 0.05)
  labels1 = tuk1$group[,3]
  labels2 = tuk2$group[,3]
  labels3 = tuk3$group[,3]
  len1 = length(labels1)
  len2 = length(labels2)
  len3 = length(labels3)
  gg = ggplot(dat, aes(x=Order, y=Property, fill= Horizon)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
    annotate("text", x = c(1:len1), y = rep(b, len1), label = labels1, 
             size = s, colour = "black") +
    annotate("text", x = c(1:len2), y = rep(c, len2), label = labels2, 
             size = s, colour = "black") +
    annotate("text", x = c(1:len3), y = rep(d, len3), label = labels3, 
             size = s, colour = "black") +
    #    annotate("text", x = c(seq(v1, v1+(len-1)*w, by = w)), y = rep(c, len), label = labels1, 
    #             size = s, colour = "black") +
    #    annotate("text", x = c(seq(v2, v2+(len-1)*w, by = w)), y = rep(b, len), label = labels2, 
    #             size = s, colour = "black") +
    #    annotate("text", x = c(seq(v3, v3+(len-1)*w, by = w)), y = rep(c, len), label = labels3, 
    #             size = s, colour = "black") +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  
  if(trans == "log"){
    gg = gg + scale_y_continuous(trans=log10_trans())
  }
  
  return(gg)
}
# For Top/Sub
setcomp2 = function(dat, unit, b, c, s, t1, t2, trans){
  # Args:
  # * dat: data only includes property, Order and Horizon by first, second and third column.(no NA)
  # * unit: the unit of property
  # * b: vertical coordinate of Tukey's comparison label.of (A Horizon)
  # * c: vertical coordinate of Tukey's comparison label.of (B Horizon)
  # * d: vertical coordinate of Tukey's comparison label.of (C Horizon)
  # * s: size of the label.
  # * t1, t2: size of axis text and axis title.
  # * trans: what transformation to the data, "log" stands for log10 transformation.
  
  input = paste0(colnames(dat)[1], unit)
  colnames(dat)[1] = "Property"
  Topdat = dat[which(dat[,3] == "Top"), ]
  Subdat = dat[which(dat[,3] == "Sub"), ]
  model1 = aov(Property~Order, data = Topdat)
  tuk1 = HSD.test(model1, "Order", group=TRUE, alpha = 0.05)
  model2 = aov(Property~Order, data = Subdat)
  tuk2 = HSD.test(model2, "Order", group=TRUE, alpha = 0.05)
  labels1 = tuk1$group[,3]
  labels2 = tuk2$group[,3]
  len1 = length(labels1)
  len2 = length(labels2)
  gg = ggplot(dat, aes(x=Order, y=Property, fill= Horizon)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
    annotate("text", x = c(1:len1), y = rep(b, len1), label = labels1, 
             size = s, colour = "black") +
    annotate("text", x = c(1:len2), y = rep(c, len2), label = labels2, 
             size = s, colour = "black") +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  
  if(trans == "log"){
    gg = gg + scale_y_continuous(trans=log10_trans())
  }
  
  return(gg)
}

#Hist over each horizon in every order(set)
orderHist = function(dat, k,x, y1,y2,y3, y4, pro, wid){
  # Args:
  # * k: which horizon 1:A, 2:B, 3:C
  # * x: horizon coordinate of mean sd and N and "Horizon"
  # * y1, y2, y3: vertical coordinate of mean sd and N
  # * y4: vertical position of "A|B|C Horizon"
  # * pro: which property? eg: 2 represents sand.
  # * wid: width of histogram.
  
  dat = na.omit(dat[,c(1,pro)])
  
  horizon = levels(dat[,1])[k]
  dat = dat[which(dat[,1] == horizon),]
  dat = droporder(dat,1)
  m1 = mean(dat[,2])
  s1 = sd(dat[,2])
  arg1 = paste("mean =", round(m1, digits=3))
  arg2 = paste("sd =", round(s1, digits=3))
  arg3 = paste("N =", dim(dat)[1])
  colnames(dat)[2] = "Property"
  if(levels(dat[,1]) == "A"){
    colors = "#CC6666"
    input1 = "Relative Frequency"
  }else{
    input1 = ""
    if(levels(dat[,1]) == "B"){
      colors = "#9999CC"
    }else{
      colors = "#66CC99"
    }
  }
  input2 = paste(horizon, "Horizon")
  plot = ggplot(dat, aes(x = Property)) + 
    geom_histogram(aes(y = ..density..), binwidth= wid, fill = colors, colour = "black") + 
    ylab(input1) +
    xlab("") +
    #scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C") +
    stat_function(fun = dnorm, colour = "red",  
                  args = list(mean = m1, sd = s1)) +
    annotate("text", x = rep(x,3), y = c(y1,y2,y3), label = c(arg1, arg2, arg3),size = 7, colour = "black") +
    annotate("text", x = x, y = y4, label = input2, colour = colors, size = 10) +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=15,face="bold"))
  return(plot)
}

#if same position
library(gridExtra)
easygrid = function(dat, a,b1,b2,b3,b4, pro, wid){
  a1 = orderHist(dat,1,a,b1,b2,b3,b4,pro,wid)
  a2 = orderHist(dat,2,a,b1,b2,b3,b4,pro,wid)
  a3 = orderHist(dat,3,a,b1,b2,b3,b4,pro, wid)
  name = colnames(dat)[pro]
  tables = table(dat[,11])
  title = paste(name, "in", names(which(tables>0)))
  par(mfrow=c(1,3))
  #grid.arrange(a1, a2, a3, ncol = 3, nrow = 1, top = title)
  grid.arrange(a1, a2, a3, ncol = 3, nrow = 1)
}

library(grid)
library(gridExtra)

# tansform dataset to have Property columns in order to plot circle
transdat = function(dat, m){
  dat = na.omit(dat[, c(1, m, 11)])
  if(colnames(dat)[2] == "Organic.Carbon"){
    dat = dat[-which(dat[,2]<0),]
  }else{
    if(colnames(dat)[2] == "Silt")
      dat = dat[-which(dat[,2]<0), ]
    }
  order = paste(dat[,3], dat[,1])
  dat = cbind(dat, order =order)
  dat = dat[, -c(1,3)]
  n = length(levels(dat[,2]))
  le = levels(dat[,2])
  num = as.numeric(n)
  for(i in 1:n){
    num[i] = mean(dat[which(dat[,2] == le[i]),1])
  }
  num = round(num, digits = 3)
  nums = data.frame(value = num, property = rep(colnames(dat)[1], n), order = le)
  return(nums)
}

# generate the circle point plot
circlePlot = function(dat, group){
  un1 = "(g kg^-1)"
  un3 = "(cmol kg^-1)"
  un4 = "(% of CEC)"
  if(group == "physical"){
    h1 = transdat(dat, 2)
    h2 = transdat(dat, 3)
    h3 = transdat(dat, 4)
    h4 = transdat(dat, 9)
    h1 = addUnit(h1, 2, "(%)")
    h2 = addUnit(h2, 2, "(%)")
    h3 = addUnit(h3, 2, "(%)")
    h4 = addUnit(h4, 2, "(% of CEC)")
    total = rbind(h1,h2,h3,h4)
  }else{
    h5 = transdat(dat, 5)
    h6 = transdat(dat, 8)
    h7 = transdat(dat, 10)
    h5 = addUnit(h5, 2, un1)
    h6 = addUnit(h6, 2, un3)
    h7 = addUnit(h7, 2, un2)
    total = rbind(h5,h6,h7)
  }
  gg = ggplot(data=total, aes(x=order, y=value, group=property, colour=property)) +
    geom_line() +
    geom_point() +
    xlab("") +
    ylab("") +
    theme(axis.text=element_text(size=8, face="bold"),axis.title=element_text(size=9,face="bold")) +
    coord_polar(theta = "x", direction=1 )
  return(gg)
}

# adding unit
addUnit = function(dat, m, unit){
  temp = as.character(dat[, m])
  temp = paste(temp, unit)
  temp = as.factor(temp)
  dat[, m] = temp
  return(dat)
}

# Drop None(Zero) Orders in a column.
droporder = function(dat, m){
  # Args:
  # * dat: the Original dataset.
  # * m: the target column you want to drop levels.
  
  worder = dat[,m]
  worder = as.character(worder)
  worder = as.factor(worder)
  dat[,m] = worder
  sign = is.na(worder)
  dat = dat[!sign,]
  return(dat)
}

###Clustering####
# Determine number of clusters
num_clusters = function(dat, m){
  # Args:
  # * m: number of clusters you want to try.
  
  hwss <- (nrow(dat)-1)*sum(apply(dat,2,var))
  for (i in 2:m) {
    hwss[i] <- sum(kmeans(dat,centers=i, nstart = 20)$withinss)
  }
  gg = plot(1:m, hwss, type="b", xlab="Number of Clusters",
            ylab="Within groups sum of squares")
  return(gg)
}

# Boxplot for data after k-means method (Including Tukey pairwise comparison)
library(scales)
ClustersComp = function(dat, unit, v, s, t1, t2, trans){
  # Args:
  # * dat: data only includes property, Order and Horizon by first, second and third column.(no NA)
  # * unit: the unit of property
  # * v: vertical coordinate of Tukey's comparison label.
  # * s: size of the label.
  # * t1, t2: size of axis text and axis title.
  # * trans: what kind of transfomation to y axis.
  
  input = paste(colnames(dat)[2], unit)
  colnames(dat)[2] = "Property"
  model = aov(Property~Clusters, data = dat)
  tuk = HSD.test(model, "Clusters", group = TRUE, alpha = .05)
  labels = tuk$group[,3]
  len = length(labels)
  tables = table(dat[,1])
  ho = as.character(dat[,1])
  ho[which(ho == "1")] = paste0("1 (n = ",tables[1], ")")
  ho[which(ho == "2")] = paste0("2 (n = ",tables[2], ")")
  ho[which(ho == "3")] = paste0("3 (n = ",tables[3], ")")
  ho[which(ho == "4")] = paste0("4 (n = ",tables[4], ")")
  ho[which(ho == "5")] = paste0("5 (n = ",tables[5], ")")
  ho[which(ho == "6")] = paste0("6 (n = ",tables[6], ")")
  ho[which(ho == "7")] = paste0("7 (n = ",tables[7], ")")
  ho = as.factor(ho)
  dat[,1] = ho
  gg = ggplot(dat, aes(x=Clusters, y=Property, fill= Clusters)) + geom_boxplot() +
    ylab(input)+
    xlab("")+
    guides(fill = F) +
    annotate("text", x = c(1:len), y = rep(v, len), label = labels, 
             size = s, colour = "black") +
    theme(axis.text=element_text(size=t1),axis.title=element_text(size=t2,face="bold"))
  
    
  if(trans == "log"){
    gg = gg + scale_y_continuous(trans=log10_trans())
  }
  
  return(gg)
}

#TukeyBoxplot2 = function(dat, a, b, a1, b1, b2, unit, s1){
  # Args:
  # * a,b: number of columns selected in dat.
  # * a1, b1, b2, s1: x axis and y axis, also the by in seq.for p-vlaue
  # * unit: unit of property
  
  
  dat = dat[, c(a,b)]
  input3 = paste(colnames(dat)[1], unit)
  colnames(dat)[1] = "Property"
  model <- aov(Property ~ Order, data = dat)
  result = TukeyHSD(model)$Order
  temp = result[,4]
  if(all(temp < 0.05)){
    gg = ggplot(dat, aes(x=Order, y=Property, fill=Order)) + geom_boxplot()+ 
      ylab(input3) +
      xlab("") +
      guides(fill = F) +
      stat_summary(fun.y= mean, colour="darkred", geom="point", shape=18, size=3, show_guide = FALSE) +
      stat_summary(fun.y= mean, colour="black", geom="text", size = 8, show_guide = FALSE, 
                   vjust=-0.7, aes( label=round(..y.., digits=1))) +
      stat_summary(fun.y=max, colour="black", geom="text", size = 8, show_guide = FALSE, 
                   vjust=-0.2, aes( label=round(..log10(y).., digits=1))) +
      stat_summary(fun.y=min, colour="black", geom="text", size = 8, show_guide = FALSE, 
                   vjust=-0.7, aes( label=round(..y.., digits=1))) +
      theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"))
  }else{
    tt = result[result[,4]>0.05,4]
    nn =  names(which(result[,4] > 0.05))
    tt = round(tt, digits = 3)
    p = NULL
    g = " = "
    for(i in 1:length(tt)){
      p[i] = paste0(nn[i], g)
      p[i] = paste0(p[i], tt[i])
    }
    len = length(p)
    gg = ggplot(dat, aes(x=Order, y=Property, fill=Order)) + geom_boxplot()+ 
      ylab(input3) +
      xlab("") +
      guides(fill = F) +
      stat_summary(fun.y= mean, colour="darkred", geom="point", shape=18, size=3, show_guide = FALSE) +
      stat_summary(fun.y= mean, colour="black", geom="text", size = 8, show_guide = FALSE, 
                   vjust=-0.7, aes( label=round(..y.., digits=1))) +
      stat_summary(fun.y=max, colour="black", geom="text", size = 8, show_guide = FALSE, 
                   vjust=-0.2, aes( label=round(..y.., digits=1))) +
      stat_summary(fun.y=min, colour="black", geom="text", size = 8, show_guide = FALSE, 
                   vjust=-0.7, aes( label=round(..y.., digits=1))) +
      theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold")) +
      annotate("text", x = rep(a1,len), y = seq(b1, b1+(len-1)*b2, by = b2), label = p, colour = "blue", size = s1)
    
  }
  return(gg)
  
}


###Classification Method Comparisons###
library("ggplot2")
classification.comparisons = function(dat, m, k, u){
  # Args:
  #  * m: mtry entry
  #  * k: setseed number
  #  * u: cols of dat - 1 (how many predictors)
  # return the 5 missclassification rate
  miss = numeric(5)
  set.seed(k)
  miss[1] = cv.rf(dat, m)
  miss[2] = cv.logit(dat)
  miss[3] = cv.tree(dat)
  miss[4] = cv.lda(dat, u)
  miss[5] = cv.svm(dat)
  miss = miss * 100
  return(miss)
  
}
# plot the comparison of classification rate
ratecomp = function(miss, wid, a1, a2, a3){
  # Args:
  # * wid: width of bar.
  # * a1, a2, a3: size of max, axis and axis title
  
  rates = data.frame(Method = c("randomForest", "LinearDiscriminantAnalysis", "Logistic", "Tree", 
                                "Support Vector Machine"), error_rate = miss)
  gg = ggplot(data=rates, aes(x=Method, y=error_rate, fill=Method)) +
    geom_bar(stat="identity", width = wid)  +
    stat_summary(fun.y=max, colour="black", geom="text", size = a1, show_guide = FALSE, 
                 vjust=-0.2, aes( label=round(..y.., digits=2))) +
    xlab("")+
    ylab("misclassification rate (%)")+
    theme(axis.text=element_text(size=a2),axis.title=element_text(size=a3,face="bold"))
  return(gg)
}
# cv on randomForest
library("randomForest")
cv.rf = function(dat, m){
  index = sample(nrow(dat))
  group = rep(1:5, length = nrow(dat))
  f_index = split(index, group)
  miss = as.numeric(5)
  set.seed(1)
  for(i in 1:5){
    datas = dat[-f_index[[i]], ]
    folds = dat[f_index[[i]], ]
    model = randomForest(Order~., data = datas, mtry = m)
    pred = predict(model, folds)
    miss[i] = 1 - sum(diag(table(pred, folds$Order)))/(dim(folds)[1])
  }
  ave_miss = mean(miss)
  return(ave_miss)
}
# cv on Logistic Regression.
library("nnet")
cv.logit = function(dat){
  index = sample(nrow(dat))
  group = rep(1:5, length = nrow(dat))
  f_index = split(index, group)
  miss = as.numeric(5)
  for(i in 1:5){
    datas = dat[-f_index[[i]], ]
    folds = dat[f_index[[i]], ]
    model = multinom(formula = Order~., data = datas)
    pred = predict(model, folds)
    miss[i] = 1 - sum(diag(table(pred, folds$Order)))/(dim(folds)[1])
  }
  ave_miss = mean(miss)
  return(ave_miss)
}

# cv on Tree Method.
library("rpart")
cv.tree = function(dat){
  index = sample(nrow(dat))
  group = rep(1:5, length = nrow(dat))
  f_index = split(index, group)
  miss = as.numeric(5)
  for(i in 1:5){
    datas = dat[-f_index[[i]], ]
    folds = dat[f_index[[i]], ]
    model = rpart(formula = Order~., data = datas)
    pred = predict(model, folds, type = "class")
    miss[i] = 1 - sum(diag(table(pred, folds$Order)))/(dim(folds)[1])
  }
  ave_miss = mean(miss)
  return(ave_miss)
}

# cv on Linear Discriminant Analysis.
library(classifly)
library(MASS)
cv.lda = function(dat, m){
  index = sample(nrow(dat))
  group = rep(1:5, length = nrow(dat))
  f_index = split(index, group)
  miss = as.numeric(5)
  for(i in 1:5){
    datas = dat[-f_index[[i]], ]
    folds = dat[f_index[[i]], ]
    x1 = as.matrix(datas[, 1:m])
    x2 = as.matrix(folds[, 1:m])
    y1 = as.numeric(datas[, m+1])
    y2 = as.numeric(folds[, m+1])
    model = lda(x1, y1)
    pred = predict(model, x2)$class
    miss[i] = 1 - sum(diag(table(pred, y2)))/(dim(folds)[1])
  }
  ave_miss = mean(miss)
  return(ave_miss)
}

#cv on support vector machine
library("e1071")
cv.svm = function(dat){
  index = sample(nrow(dat))
  group = rep(1:5, length = nrow(dat))
  f_index = split(index, group)
  miss = as.numeric(5)
  for(i in 1:5){
    datas = dat[-f_index[[i]], ]
    folds = dat[f_index[[i]], ]
    model = svm(Order~., data=datas, kernel="radial", cost=10, gamma=1)
    pred = predict(model, folds)
    miss[i] = 1 - sum(diag(table(pred, folds$Order)))/(dim(folds)[1])
  }
  ave_miss = mean(miss)
  return(ave_miss)
}

# Variables Importance Plot on randomForest.
rf.Imp = function(dat, m, k){
  # Args:
  # * dat: data that you want to do the classification.
  # * m: randomForest mtry entry.
  # * k: set.seed number.
  
  set.seed(k)
  model = randomForest(Order~.,data=dat, importance = T, mtry=m)
  #importance(model)
  #varImpPlot(model)
  
  imp <- importance(model, type=1)
  featureImportance <- data.frame(Feature=row.names(imp), Importance=imp[,1])
  
  p <- ggplot(featureImportance, aes(x=reorder(Feature, Importance), y=Importance)) +
    geom_bar(stat="identity", fill="#53cfff") +
    coord_flip() + 
    theme_light(base_size=20) + 
    xlab("") +
    ylab("") +
    theme(plot.title=element_text(size=18))
  
  return(p)
}