package ch.irb.ManageFastaFiles;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;


/**
 * Copyright 2017 - Mathilde Foglierini Perez
 * This code is distributed open source under the terms of the GNU Free Documention License.
 * <p>
 * This class will create a fasta file. The input will be the path of the file, then  a hashMap
 * containing all the IDs with their sequences
 */
public class FastaFileMaker {

    private String filePath;
    private HashMap<String, String> fastaIdToSequence = new HashMap<String, String>();

    public FastaFileMaker(String filePath, HashMap<String, String> fastaIdToSequence) {
        setFilePath(filePath);
        setFastaIdTosequence(fastaIdToSequence);
        try {
            createFastaFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void createFastaFile() throws IOException {
        try {
            FileWriter fstream = new FileWriter(filePath);
            BufferedWriter out = new BufferedWriter(fstream);
            for (Map.Entry<String, String> entry : fastaIdToSequence.entrySet()) {
                String fastaId = entry.getKey();
                String sequence = entry.getValue();
                out.write(">" + fastaId
                        + System.getProperty("line.separator") + sequence
                        + System.getProperty("line.separator"));
            }
            // Close the output stream
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public String getFilePath() {
        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }

    public HashMap<String, String> getFastaIdTosequence() {
        return fastaIdToSequence;
    }

    public void setFastaIdTosequence(HashMap<String, String> fastaIdTosequence) {
        this.fastaIdToSequence = fastaIdTosequence;
    }

}
