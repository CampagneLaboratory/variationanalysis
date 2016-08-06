# Use this script to look at features that are very different between two groups of predictions (i.e., correct, wrong).
require("data.table")
require("dtplyr")
require("tidyr")

combined<-fread(header=TRUE,"/Users/fac2003/R_RESULTS/variation/combined.tsv",colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","string"))

small <-sample_n(combined, size=10000)

tidy <- combined %>%
    group_by(GROUP) %>%
    gather( na.rm=TRUE,  feature, value, -id, -sumCounts, -GROUP) %>%
    select(id, sumCounts, feature, value, GROUP)
head(tidy, n=10)

full<- tidy %>% group_by(GROUP, sumCounts, feature) %>%
summarize( n=n(), mean=mean(value,na.rm = TRUE), sd=sd(value, na.rm = TRUE)) %>%
filter(n>=50)
head(tidy, n=10)

correct <- full %>% filter(GROUP=="Correct")
wrong <- full %>% filter(GROUP=="Wrong")
joined<-inner_join(correct,wrong, by=c("feature" = "feature"))
head(joined,n=10)

sorted <-  joined%>% group_by(feature) %>%
mutate(meanCorrect=mean.x,
meanWrong=mean.y,
diff=abs(meanWrong-meanCorrect)) %>%
select( diff,feature) %>%
group_by(feature) %>%
summarize(maxDiff=max(diff)) %>%
arrange(desc(maxDiff))

View(head(sorted,n=10))
