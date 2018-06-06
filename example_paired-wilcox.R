## 20180524
## http://www.sthda.com/english/wiki/paired-samples-t-test-in-r


# Weight of the mice before treatment
before <-c(200.1, 190.9, 192.7, 213, 241.4, 196.9, 172.2, 185.5, 205.2, 193.7)
# Weight of the mice after treatment
after <-c(392.9, 393.2, 345.1, 393, 434, 427.9, 422, 383.9, 392.3, 352.2)
# Create a data frame
my_data <- data.frame( 
                group = rep(c("before", "after"), each = 10),
                weight = c(before,  after)
                )

libr("dplyr")
group_by(my_data, group) %>%
  summarise(
    count = n(),
    median = median(weight, na.rm = TRUE),
    IQR = IQR(weight, na.rm = TRUE)
  )

# devtools::install_github("kassambara/ggpubr")

libr("ggpubr")
ggboxplot(my_data, x = "group", y = "weight", 
          color = "group", palette = c("#00AFBB", "#E7B800"),
          order = c("before", "after"),
          ylab = "Weight", xlab = "Groups")

# Subset weight data before treatment
before <- subset(my_data,  group == "before", weight,
                 drop = TRUE)
# subset weight data after treatment
after <- subset(my_data,  group == "after", weight,
                drop = TRUE)
# Plot paired data
libr(PairedData)
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()

res <- t.test(before, after, paired = TRUE)
res <- t.test(before, after, paired = TRUE)

# test whether the average weight before treatment is less than the average weight after treatment
res <- t.test(weight ~ group, data = my_data, paired = TRUE,
       alternative = "less")

# test whether the average weight before treatment is greater than the average weight after treatment
res <- t.test(weight ~ group, data = my_data, paired = TRUE,
       alternative = "greater")

res$p.value
res$estimate
res$conf.int












## 2

res <- wilcox.test(before, after, paired = TRUE)

# Compute t-test
res <- wilcox.test(weight ~ group, data = my_data, paired = TRUE)

# print only the p-value
res$p.value

# test whether the median weight before treatment is less than the median weight after treatment
res <- wilcox.test(weight ~ group, data = my_data, paired = TRUE,
            alternative = "less")

# test whether the median weight before treatment is greater than the median weight after treatment
res <- wilcox.test(weight ~ group, data = my_data, paired = TRUE,
            alternative = "greater")