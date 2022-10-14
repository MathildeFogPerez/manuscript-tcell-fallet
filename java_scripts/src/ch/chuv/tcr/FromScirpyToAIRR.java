package ch.chuv.tcr;

import ch.irb.utils.Consts;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

/*
Copyright 2022 - Mathilde Foglierini Perez
This class will take in input the airr file produced by scirpy. It will clean the data by removing non producing TCR,
add the full VDJ nucleotide and protein sequence for each clonotype using the AIRR file produced by changeO
and split the data in X files where X is the number of donors. Those files can be used as input for immunarch for example.
It will also generate one single AIRR file.
 */
public class FromScirpyToAIRR {

    //TODO set here the right path
    private File folder = new File("/MYPATH/manuscript-tcell-fallet/analysis_out_10X/files");
    private File cellranger_folder = new File("/MYPATH/manuscript-tcell-fallet/cellranger_out");
    private File scirpyAirrFile = new File(folder.getPath()+ Consts.fs+"airrfile_data.tsv");
    private static ArrayList<String> donors = new ArrayList<>();
    private  HashMap<String,BufferedWriter> donorToBffWriter = new HashMap<>();
    private HashMap<String,String> sequenceIdToProtSeq = new HashMap<>();
    private HashMap<String,String> sequenceIdToNucSeq = new HashMap<>();

    public static void main(String[] args) throws IOException {
        //TODO here put the right sample names
        donors.add("BP_01"); donors.add("BP_02");donors.add("BP_03");
        new FromScirpyToAIRR();
    }


    public FromScirpyToAIRR() throws IOException {
        for (String donor: donors){
            parseAIRRfile(donor);
        }
        //set the files
        BufferedWriter metadataout = new BufferedWriter(new FileWriter(new File(folder.getPath()+Consts.fs+"metadata.txt")));
        metadataout.write("Sample"+Consts.ls);
        for (String donor: donors){
            File outFile = new File(folder.getPath()+Consts.fs+donor+".tsv");
            BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
            donorToBffWriter.put(donor,out);
            metadataout.write(donor+Consts.ls);
        }
        metadataout.close();
        //Write in a single file (used for search of Anti-spike TCR spe)
        BufferedWriter singleOut = new BufferedWriter(new FileWriter(new File(folder.getPath()+Consts.fs
                +scirpyAirrFile.getName().replace(".tsv","_cleaned.tsv"))));
        BufferedReader fileReader = new BufferedReader(new FileReader(scirpyAirrFile));
        String line="";
        int index=0;
        int seqIndex=1;
        HashMap<String,Integer>headerToIndex= new HashMap<>();
        while ((line= fileReader.readLine())!=null){
            String[] cells = line.split("\t");
            index++;
            if (index==1){
                int i=0;
                for (String cell:cells){
                    headerToIndex.put(cell,i);
                    i++;
                }
                singleOut.write(line);
                singleOut.write("\tsequence_aa"+Consts.ls);
                for (String donor: donors) {
                    donorToBffWriter.get(donor).write(line + Consts.ls);
                }
            }
            else{
                String isProductive= cells[3];
                if (isProductive.equals("F")){
                    continue;
                }
                String donor = cells[headerToIndex.get("sample")];
                BufferedWriter out = donorToBffWriter.get(donor);
                String seqName=donor;
                if (cells[4].contains("TRA")){
                    seqName+="_TRA";
                }else{
                    seqName+="_TRB";
                }
                seqName+="_"+cells[headerToIndex.get("cell_id")];
                out.write(seqName+"\t\tT");
                singleOut.write(seqName+"\t\tT");
                //System.out.println("Checking "+seqName);
                if (seqName.contains("-1-1")){
                    seqName= seqName.substring(0,seqName.length()-2);
                }
                for (int i=3;i<cells.length;i++){
                    if (headerToIndex.get("sequence_alignment") ==i){
                        out.write("\t"+sequenceIdToNucSeq.get(seqName));
                        singleOut.write("\t"+sequenceIdToNucSeq.get(seqName));
                    }
                    else{
                        out.write("\t"+cells[i]);
                        singleOut.write("\t"+cells[i]);
                    }
                }
                String seq= sequenceIdToProtSeq.get(seqName);
                //check we have the right sequence!
                if (!seq.contains(cells[headerToIndex.get("junction_aa")])){
                    System.out.println("GO OUT, problem with "+seqName+" and "+seq);
                    System.exit(-1);
                }
                singleOut.write("\t"+seq+Consts.ls);
                out.write(Consts.ls);
                seqIndex++;
            }
        }
        fileReader.close();
        singleOut.close();
        //close all files
        for (String donor: donors) {
            donorToBffWriter.get(donor).close();
        }
        seqIndex--;
        index--;
        System.out.println("Number of sequences parsed = "+index);
        System.out.println("Number of cleaned sequences written = "+seqIndex);
    }
    private void parseAIRRfile(String donor) throws IOException {
        File airrfile = new File(cellranger_folder.getPath()+donor+"_airr_rearrangement.tsv");
        BufferedReader fileReader = new BufferedReader(new FileReader(airrfile));
        String line="";
        int index=0;
        HashMap<String,Integer>headerToIndex= new HashMap<>();
        while ((line= fileReader.readLine())!=null) {
            String[] cells = line.split("\t");
            index++;
            if (index == 1) {
                int i=0;
                for (String cell:cells){
                    headerToIndex.put(cell,i);
                    i++;
                }
            }else{
                String seqId=donor+"_"+cells[headerToIndex.get("j_call")].substring(0,3)+"_"+cells[headerToIndex.get("cell_id")];
                //System.out.println(seqId);
                String seq= cells[headerToIndex.get("sequence_aa")];
                sequenceIdToProtSeq.put(seqId,seq);
                String dnaSeq= cells[headerToIndex.get("sequence_alignment")];
                sequenceIdToNucSeq.put(seqId,dnaSeq);
            }
        }
        fileReader.close();
    }
}
