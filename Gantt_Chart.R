# gantt chart
library(plotrix)
library(lubridate)

# set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create a list of items for the Gantt Chart
tasks=list(
  # task labels
  labels=c("Literature Review", "Presentation Preparation", "Dataset Acquisition", "Subset Identification", 
           "Ligand-Receptor Pairing", "Gene Expression Analysis", "Functional Interpretation", "Spatial Integration", "Draft Write-up", "Final Write-up"),
  # start dates of each task
  starts=ymd("2026-04-17", "2026-05-16", "2026-04-17", "2026-05-16", "2026-06-05", "2026-06-19", "2026-07-03", "2026-07-10", "2026-06-05", "2026-08-10"),
  # end dates of each task
  ends=ymd("2026-05-15", "2026-05-26", "2026-05-01", "2026-06-05", "2026-06-19","2026-07-03", "2026-07-10", "2026-07-26", "2026-08-07", "2026-08-31"),
  # "priorities"- will be used as colours for task type
  priorities=c(4, 4, 1, 1, 2, 2, 3, 3, 4, 4)
)

# set the x axis labels
vgridlab=c("April","May", "June", "July", "August", "September")
# and where they will be
vgridpos=ymd("2026-04-01", "2026-05-01", "2026-06-01", "2026-07-01", "2026-08-01", "2026-09-01")

# save pdf of the Gantt chart
pdf("gantt_chart.pdf", width = 10, height = 6)

# generate gantt chart
gantt.chart(tasks,vgridpos = vgridpos, 
            vgridlab = vgridlab, 
            # remove priority legend
            priority.legend = FALSE,
            hgrid= TRUE,
            # empty title
            main = "",
            border.col = "black",
            # "priority colours"
            taskcolors = c("forestgreen", "orange" , "brown1", "skyblue"),
            # add extra room on each side
            xlim = ymd(c("2026-04-10", "2026-09-07"))
            )

# add a legend as for task types instead of priorities
legend(x = as.Date("2026-07-22"),
       y = 10.1, 
       bg = "white",
       legend = c("Data Preparation", "Interaction Analysis", "Biological interpretation", "Report Writing"),
       fill = c("forestgreen", "orange" , "brown1", "skyblue"))

dev.off()
