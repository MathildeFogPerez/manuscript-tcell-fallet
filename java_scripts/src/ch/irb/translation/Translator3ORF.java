package ch.irb.translation;

import java.util.HashMap;




/**
 * @author Mathilde This class take the DNA sequence and translates it into a protein sequence
 */
public class Translator3ORF {
	private String dnaSequence;
	private String protSequence = "";
	private String partOfSequence;
	private boolean checkPieceOfSeq=false;
	@SuppressWarnings("unused")
	private boolean isDNA;
	private HashMap<String, String> codons = new HashMap<String, String>();

	public Translator3ORF(String dnaSequence, boolean isDNA) {
		this(dnaSequence,isDNA,null);
	}

	public Translator3ORF(String dnaSequence, boolean isDNA, String substring) {
		partOfSequence= substring;
		if (substring != null){
			checkPieceOfSeq = true;
		}
		this.isDNA = isDNA;
		dnaSequence = dnaSequence.toUpperCase();
		this.dnaSequence = dnaSequence;
		if (!isDNA){
			protSequence = dnaSequence ;
		}else{
			setCodonsTable();
			try {
				//translate();
				String prot1 =translate(0);
				String prot2 =translate(1);
				String prot3 =translate(2);
				protSequence=prot1;
				//logger.warn("prot1: "+prot1);
				
				//bug fixed the 20.07.20
				//remove  && !prot2.contains("-")
				if (!prot2.contains("*") && (prot1.contains("*"))){
					protSequence=prot2;
					//logger.warn("prot2: "+prot2);
				}
				//bug fixed the 20.07.20
				//remove  && !prot3.contains("-")
				else if (!prot3.contains("*") && (prot1.contains("*"))){
					protSequence=prot3;
					//logger.warn("prot3: "+prot3);
				}
				
				//if all 3 proteins have stop codon, we count the number of stop codons 
				//and take the protein where there are the less
				else if (prot1.contains("*") && prot2.contains("*") && prot3.contains("*")){
					int stopNumbProt1 = 0;
					for (char aa: prot1.toCharArray()){
						if (aa == '*'){
							stopNumbProt1++;
						}
					}
					int stopNumbProt2 = 0;
					for (char aa: prot2.toCharArray()){
						if (aa == '*'){
							stopNumbProt2++;
						}
					}
					int stopNumbProt3 = 0;
					for (char aa: prot3.toCharArray()){
						if (aa == '*'){
							stopNumbProt3++;
						}
					}
					//by default it will be prot 1
					if (stopNumbProt2 < stopNumbProt1 && stopNumbProt2 < stopNumbProt3){
						protSequence = prot2;
					}
					else if (stopNumbProt3 < stopNumbProt1 && stopNumbProt3 < stopNumbProt2){
						protSequence = prot3;
					}
				}
				//special case when we know a part of the sequence (i.e. FWR4 for Ig)
				if(checkPieceOfSeq){
					if (prot1.contains(partOfSequence)){
						protSequence=prot1;
					}
					else if (prot2.contains(partOfSequence)){
						protSequence=prot2;
					}
					else if (prot3.contains(partOfSequence)){
						protSequence=prot3;
					}
				}
				
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
	}

	private String translate(int start) throws Exception {
		// logger.warn("Sequence to translate is "+dnaSequence);
		String prot = "";
		char[] nuc = new char[dnaSequence.length()];
		nuc = dnaSequence.toCharArray();
		String codon = "";
		int nucleotide = 0;
		for (int i = start; i < nuc.length; i++) {
			nucleotide++;
			if (nucleotide == 1 || nucleotide == 2 || nucleotide == 3) {
				codon = codon.concat(String.valueOf(nuc[i]));
			}
			if (nucleotide == 3) {
				//logger.warn("CODON is " + codon);
				// here we set the prot to null if there is indels that doesnt allow the translation
				// TODO here deal with the 3 ORF
				String aa;
				if (!codons.containsKey(codon)) { //case of N M K R S Y W nuc from sequencing data
					aa = "X";
				} else {
					aa = codons.get(codon);
				}
				prot = prot.concat(aa);
				// logger.warn("We have the codon "+codon+ " that codes for "+aa);
				// logger.debug("The protein sequence is "+protSequence);
				nucleotide = 0;
				codon = "";

			}
		}
		return prot;
	}
	

	private void setCodonsTable() {
		codons.put("ATT", "I");
		codons.put("ATC", "I");
		codons.put("ATA", "I");
		codons.put("CTT", "L");
		codons.put("CTC", "L");
		codons.put("CTA", "L");
		codons.put("CTG", "L");
		codons.put("TTA", "L");
		codons.put("TTG", "L");
		codons.put("GTT", "V");
		codons.put("GTC", "V");
		codons.put("GTA", "V");
		codons.put("GTG", "V");
		codons.put("TTT", "F");
		codons.put("TTC", "F");
		codons.put("ATG", "M");
		codons.put("TGT", "C");
		codons.put("TGC", "C");
		codons.put("GCT", "A");
		codons.put("GCC", "A");
		codons.put("GCA", "A");
		codons.put("GCG", "A");
		codons.put("GGT", "G");
		codons.put("GGC", "G");
		codons.put("GGA", "G");
		codons.put("GGG", "G");
		codons.put("CCT", "P");
		codons.put("CCC", "P");
		codons.put("CCA", "P");
		codons.put("CCG", "P");
		codons.put("ACT", "T");
		codons.put("ACA", "T");
		codons.put("ACG", "T");
		codons.put("ACC", "T");
		codons.put("TCT", "S");
		codons.put("TCC", "S");
		codons.put("TCA", "S");
		codons.put("TCG", "S");
		codons.put("AGT", "S");
		codons.put("AGC", "S");
		codons.put("TAT", "Y");
		codons.put("TAC", "Y");
		codons.put("TGG", "W");
		codons.put("CAA", "Q");
		codons.put("CAG", "Q");
		codons.put("AAT", "N");
		codons.put("AAC", "N");
		codons.put("CAT", "H");
		codons.put("CAC", "H");
		codons.put("GAA", "E");
		codons.put("GAG", "E");
		codons.put("GAT", "D");
		codons.put("GAC", "D");
		codons.put("AAA", "K");
		codons.put("AAG", "K");
		codons.put("CGT", "R");
		codons.put("CGC", "R");
		codons.put("CGA", "R");
		codons.put("CGG", "R");
		codons.put("AGA", "R");
		codons.put("AGG", "R");
		codons.put("TAA", "*"); // X=stop codon!
		codons.put("TAG", "*"); // X=stop codon!
		codons.put("TGA", "*"); // X=stop codon!
		codons.put("---", "-"); // DELETION
	}

	public String getProteinSequence() {
		return protSequence;
	}

}
