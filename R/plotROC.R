#
#  Copyright (C) 2017 - Garvan Institute of Medical Research
#

plotROC <- function(seqs,
                    score,
                    group,
                    label,
                    refGroup,
                    title=NULL,
                    legTitle=NULL)
{
    data <- data.frame(row.names = seqs,
                       label = label,
                       score = score,
                       group = group)
    
    data <- data[!is.na(data$score),]
    data <- data[!is.na(data$group),]
    data <- data[order(row.names(data)),]
    
    data$label <- revalue(data$label, c('FP'='0', 'TP'='1'))
    
    data <- data[data$label=='1' | data$label=='0',]
    ROC  <- NULL
    AUC  <- NULL
    
    groups <- sort(data$group)
    uniqs  <- unique(groups)
    
    if (!is.null(refGroup))
    {
        uniqs <- uniqs[uniqs != refGroup & !(uniqs %in% refGroup)]
    }
    
    stopifnot(length(uniqs) > 0)
    
    # For each non-reference group...
    for (i in c(1:length(uniqs)))
    {
        # Query group
        group <- uniqs[[i]]
        
        if (is.null(refGroup))
        {
            t <- data[data$group==group,]
        }
        else
        {
            t <- data[data$group==group | data$group==refGroup,]
        }
        
        # No FP or TP?
        if (length(unique(t$label)) == 1)
        {
            # No TP... Add a TP...
            if (unique(t$label) == '0')
            {
                t <- rbind(t, data.frame(label='1', score=0, group=group))
            }
            
            # No FP... Add a FP...
            else
            {
                t <- rbind(t, data.frame(label='0', score=0, group=group))
            }
        }
        
        t <- t[with(t, order(score)),]
        
        label <- ifelse(t$label == '1', 2, 1)
        preds <- prediction(t$score, label, label.ordering=c(1,2))
        perf  <- performance(preds, 'tpr', 'fpr')
        auc   <- performance(preds, 'auc')
        
        AUC <- rbind(AUC, data.frame(Group = group,
                                     AUC = round(unlist(auc@y.values), 4)))
        ROC <- rbind(ROC, data.frame(FPR = unlist(perf@x.values),
                                     TPR = unlist(perf@y.values), group=group))
    }
    
    ROC$group = as.factor(ROC$group)
    
    p <- ggplot(data=ROC, aes_string(x='FPR', y='TPR'))             + 
           geom_abline(intercept=0, slope=1, linetype=2)            +
           geom_path(size=1, aes_string(colour='group'), alpha=0.5) +
           labs(colour=legTitle)                                    +
           theme_bw()                                               +
           theme(plot.title = element_text(hjust = 0.5))
    
    if (!is.null(title)) { p <- p + ggtitle(title) }
    
    rownames(AUC) <- NULL
    colnames(AUC) <- c('', 'AUC')
    AUC <- AUC[order(-AUC$AUC),]
    
    print(kable(AUC, row.names=FALSE))
    print(.transformPlot(p))

    list(AUC=AUC)
}