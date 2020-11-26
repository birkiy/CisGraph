


# def linkFrom(wildcards):
#         if wildcards.raw.find("control") > 0:
#                 raw = wildcards.raw.rsplit(".", 1)[0]
#                 SRR = controlDF.loc[
#                         (controlDF["Raw"] == raw) &
#                         (controlDF["Run"] == wildcards.run),
#                         "SRR"].to_list()
#         elif wildcards.raw.find("control") == -1:
#                 SRR = sampleDF.loc[
#                         (sampleDF["Raw"] == wildcards.raw) &
#                         (sampleDF["Run"] == wildcards.run),
#                         "SRR"].to_list()
#         if wildcards.run == "U" or wildcards.run == "R1":
#                 return f"raw/{SRR[0]}_1.fastq.gz"
#         elif  wildcards.run == "R2":
#                 return f"raw/{SRR[0]}_2.fastq.gz"




def linkFrom(wildcards):
        """
        Upstream is concat
        """
        if wildcards.raw.find("control") > 0:
                raw = wildcards.raw.rsplit(".", 1)[0]
                SRX = controlDF.loc[
                        (controlDF["Raw"] == raw) &
                        (controlDF["Run"] == wildcards.run),
                        "SRX"].to_list()
        elif wildcards.raw.find("control") == -1:
                SRX = sampleDF.loc[
                        (sampleDF["Raw"] == wildcards.raw) &
                        (sampleDF["Run"] == wildcards.run),
                        "SRX"].to_list()
        return f"raw/{SRX[0]}.{{run}}.fastq.gz"




rule Links:
        input:
                linkFrom
        output:
                linkTo="links/{raw}.{run}.fastq.gz"
        message:
                "Executing Links rule from {input} to {output.linkTo}"
        shell:
                """
                ln -s ../{input} {output.linkTo}
                """



#
# def linkFromControl(wildcards):
#         control = controlDF.loc[
#                 (controlDF["Raw"] == wildcards.control),
#                 "Raw"].to_list()[0]
#         if wildcards.run == "U" or wildcards.run == "R1":
#                 return f"raw/{control}_1.fastq.gz"
#         elif  wildcards.run == "R2":
#                 return f"raw/{control}_2.fastq.gz"
#
#
# rule LinksControl:
#         input:
#                 linkFromControl
#         output:
#                 linkTo="links/{control}.control.{run}.fastq.gz"
#         message:
#                 "Executing Links rule from {input} to {output.linkTo}"
#         shell:
#                 """
#                 ln -s ../{input} {output.linkTo}
#                 """
