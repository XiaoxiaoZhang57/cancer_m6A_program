library(ggplot2)

# 示例数据
data <- data.frame(
  tissue = c("baixuebing", "bone", "fei", "gan", "guan", "mianyi", "nianyeJiangye", "pifu", "ruanzuzhi", "shen", "xianti", "xingbie", "zhongshu"),
  site = c(232367, 33675, 491853, 414982, 1013454, 60822, 477, 1579843, 41242, 159464, 907210, 699290, 208616)
)

# 配色方案
colors <- c("#947A6D", "#D7B98E", "#a6bce3", "#808F7C", "#019a99", "#0077b0", "#ffba4d", "#282152", "#caa59a", "#fdbf6f", "#69B0AC", "#EEA9A9", "#854836")

# 计算百分比
data$percentage <- data$site / sum(data$site) * 100

# 绘制饼图
p <- ggplot(data, aes(x = "", y = site, fill = tissue)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = colors) +
  labs(title = "Tissue Distribution", fill = "Tissue") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_stack(vjust = 0.5),
            size = 3)  # 调整文本大小

print(p)
ggsave(p, file='tissue site pie.pdf', width=5, height=4)
