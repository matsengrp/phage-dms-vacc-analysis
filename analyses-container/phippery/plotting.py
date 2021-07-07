

def plot_sample_binding(
        ds, 
        data_table="cpm_pre_tech", 
        show=True, 
        ymax=20000, 
        ymin=-20000, 
        saveto="group_means.png"
):
    
    y = ymax - (0.1*ymax)
    N = len(set(dtpa_df["sample_id"]))
    dtpa_df["Loc"] = dtpa_df["Loc"].astype(int)
    f = (ggplot()

        + geom_line(
                aes(x="Loc", y=data_table,group="sample_id", color="sample_id"), 
                show_legend=True, alpha=0.8, size=0.8,# color="black",
            data=dtpa_df)

        + annotate("rect",xmin=29, xmax=294, ymin=ymin, 
            ymax=ymax,alpha=0.15,fill="blue", color="black", size=0.3)
        + annotate("rect",xmin=332, xmax=523, ymin=ymin, 
            ymax=ymax,alpha=0.15,fill="aqua", color="black", size=0.3)
        + annotate("rect",xmin=530, xmax=681, ymin=ymin, 
            ymax=ymax,alpha=0.15,fill="lawngreen", color="black", size=0.3)
        + geom_vline(xintercept=685, size=0.5, color="red", linetype="dashed")
        + geom_vline(xintercept=816, size=0.5, color="red", linetype="dashed")
        + annotate("rect",xmin=816, xmax=833, ymin=ymin, 
            ymax=ymax,alpha=0.15,fill="plum", color="black", size=0.3)
        + annotate("rect",xmin=910, xmax=987, ymin=ymin, 
            ymax=ymax,alpha=0.15,fill="khaki", color="black", size=0.3)
        + annotate("rect",xmin=1162, xmax=1204, ymin=ymin, 
            ymax=ymax,alpha=0.15,fill="red", color="black", size=0.3)
        + scale_y_continuous(limits=[ymin, ymax])

        + annotate("text",x = 70, y=y, label= "NTD")
        + annotate("text",x = 372, y=y, label= "RBD")
        + annotate("text",x = 570, y=y, label= "CTD")
        + annotate("text",x = 685, y=y, label= "S1/S2")
        + annotate("text",x = 850, y=y, label= "FP")
        + annotate("text",x = 950, y=y, label= "HR1")
        + annotate("text", x = 1202, y=y, label= "HR2")

        + theme_classic()
        + theme(
            figure_size=[9, 4],
            panel_grid_major_x = element_line(size=0.5, color="black"),
            panel_grid_major_y = element_line(size=3),
            panel_grid_minor_x = element_line(size=3),
            text=element_text(size=12),
        )

        + labs(title=f"",
             y=f"Mean {data_table} (N={N})",
             x="Loc")
    )

    if saveto is not None: f.save(saveto)
    if show: f.draw()
