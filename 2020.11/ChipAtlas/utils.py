



def getLib(wildcards):
        """
        This function is used by MappingBowtie2's getFastqs and Bigwig rules' getParamsBamCom, getParamsBamCov helper functions.
        if wildcards has length 3 it came from bigwig rule otherwise mapping.
        """
        # print(wildcards, "getlib wildcards")
        # print(type(wildcards), len(wildcards))
        if wildcards.raw.find("control") > 0:
                # print("Enter getlib Control")
                raw = wildcards.raw.rsplit(".", 1)[0]
                if len(wildcards) == 3:
                        # print("Enter getlib Control stage")
                        if wildcards.stage == "merged":
                                # print("Enter getlib Control stage merged")
                                lib = controlDF.loc[controlDF["ControlName"] == raw, "Library"].to_list()[0]
                        else:
                                # print("Enter getlib Control stage merged else")
                                lib = controlDF.loc[controlDF["Raw"] == raw, "Library"].to_list()[0]
                else:
                        # print("Enter getlib Control stage else")
                        lib = controlDF.loc[controlDF["Raw"] == raw, "Library"].to_list()[0]

        else:
                # print("Enter getlib Samples")
                if len(wildcards) == 3:
                        # print("Enter getlib Samples stage")
                        if wildcards.stage == "merged":
                                # print("Enter getlib Samples stage merged")
                                lib = sampleDF.loc[sampleDF["SampleName"] == wildcards.raw, "Library"].to_list()[0]
                        else:
                                # print("Enter getlib Samples stage merged else")
                                lib = sampleDF.loc[sampleDF["Raw"] == wildcards.raw, "Library"].to_list()[0]
                else:
                        # print("Enter getlib Samples stage else")
                        lib = sampleDF.loc[sampleDF["Raw"] == wildcards.raw, "Library"].to_list()[0]
        return lib


################################################################################

def getFastqs(wildcards):
        """
        Used by:
                -> MappingBowtie2
        """
        outputD = dict()
        lib = getLib(wildcards)
        if lib == "Single":
                outputD["U"] = f"links/{wildcards.raw}.U.fastq.gz"
        elif lib == "Paired":
                outputD["R1"] = f"links/{wildcards.raw}.R1.fastq.gz"
                outputD["R2"] = f"links/{wildcards.raw}.R2.fastq.gz"
        return outputD



################################################################################


def getRepsToPool(wildcards):
        """
        Used by:
                -> Pool
        """
        if wildcards.sampleName.find("control") > 0:
                sampleName = wildcards.sampleName.rsplit(".", 1)[0]
                #print("ibne", sampleName, wildcards)
                reps = list(controlDF.loc[controlDF["ControlName"] == sampleName, "Replicate"])
                return expand("results/mapping/final/{sampleName}.{rep}.control.bam", rep=reps, sampleName=[sampleName])
        else:
                reps = list(sampleDF.loc[sampleDF["SampleName"] == wildcards.sampleName, "Replicate"])
                return expand("results/mapping/final/{{sampleName}}.{rep}.bam", rep=reps)


################################################################################


def getControl(wildcards):
        """
        Note that this function only calls one merged control bam.
        control=GSE114063, so it should be called from merged stage.
        Sorry for hard coding.
        """
        if wildcards.stage == "merged":
                control = sampleDF.loc[sampleDF["SampleName"] == wildcards.raw, "ControlName"].to_list()[0]
        else:
                control = sampleDF.loc[sampleDF["Raw"] == wildcards.raw, "ControlName"].to_list()[0]
        return f"results/mapping/merged/{control}.control.bam"


################################################################################


def linkFrom(wildcards):
        if wildcards.raw.find("control") > 0:
                raw = wildcards.raw.rsplit(".", 1)[0]
                SRR = controlDF.loc[
                        (controlDF["Raw"] == raw) &
                        (controlDF["Run"] == wildcards.run),
                        "SRR"].to_list()
        elif wildcards.raw.find("control") == -1:
                SRR = sampleDF.loc[
                        (sampleDF["Raw"] == wildcards.raw) &
                        (sampleDF["Run"] == wildcards.run),
                        "SRR"].to_list()
        if wildcards.run == "U" or wildcards.run == "R1":
                return f"raw/{SRR[0]}_1.fastq.gz"
        elif  wildcards.run == "R2":
                return f"raw/{SRR[0]}_2.fastq.gz"



################################################################################


def getAll(wildcards):
        if wildcards.samples == "merged":
                samples = sorted(set((sampleDF["SampleName"]).to_list() + (controlDF["ControlName"] + ".control").to_list()))
        elif wildcards.samples == "final":
                samples = sorted(set(list(sampleDF["Raw"]) + list(controlDF["Raw"] + ".control")))
        return expand("results/mapping/{{samples}}/{sample}.bam", sample=samples)



################################################################################
#
#
# def getParamsBamCom(wildcards):
#         """
#         Two functions are reduntand. Find a way to use only one function that can
#         be used both by bamCoverage and bamCompare.
#         """
#         lib = getLib(wildcards)
#         if lib == "Single":
#                 extend = "--extendReads 150"
#         elif lib == "Paired":
#                 extend = "--extendReads"
#         defaultParams = f"--samFlagInclude 64 {extend} --centerReads --blackListFileName {blacklistFile}"
#         if wildcards.type == "unique":
#                 return f"--ignoreDuplicates {defaultParams}"
#         elif wildcards.type == "raw":
#                 return ""
#         elif wildcards.type == "SES":
#                 return "--scaleFactorsMethod SES -l 750"
#         elif wildcards.type in ["RPKM", "BPM", "CPM"]:
#                 return f"--scaleFactorsMethod None --normalizeUsing {wildcards.type}"
#
#
# ################################################################################
#
#
# def getParamsBamCov(wildcards):
#         lib = getLib(wildcards)
#         if lib == "Single":
#                 extend = "--extendReads 150"
#         elif lib == "Paired":
#                 extend = "--extendReads"
#         defaultParams = f"--samFlagInclude 64 {extend} --centerReads --blackListFileName {blacklistFile}"
#         if wildcards.type == "unique":
#                 return f"--ignoreDuplicates {defaultParams}"
#         elif wildcards.type == "raw":
#                 return ""
#         elif wildcards.type in ["RPKM", "BPM", "CPM"]:
#                 return f"--normalizeUsing {wildcards.type}"
