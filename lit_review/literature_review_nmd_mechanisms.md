# NMD Mechanisms: Literature Review in Context of ORF Model Findings

> **About this expanded version.** Every specific mechanistic claim below is followed by a verbatim quotation from the primary literature, shown as an indented blockquote. Quotations are drawn from Results or Discussion sections of the papers wherever possible (not abstracts). PDFs for all 39 cited references are available in the adjacent `papers/` folder. Four citation errors in the original bibliography have been corrected — see the **Citation corrections** section at the end.

---

## Summary

Our deep learning model for NMD prediction reveals three principal findings that connect to distinct bodies of NMD literature: (1) PTC + downstream EJC is the dominant predictor of NMD visibility, (2) approximately 15% of NMD isoforms appear driven by 5'UTR features consistent with aberrant translation initiation, and (3) 3'UTR length does not strongly predict NMD susceptibility in our system. We additionally propose that ribosomal ORF selection may be stochastic, with NMD susceptibility reflecting the interaction between start codon strength and the NMD potential of the engaged reading frame. Here we review how these findings relate to the current state of knowledge.

---

## 1. PTC + EJC as the dominant NMD mechanism

The model's reliance on downstream exon-junction complexes as the primary predictor of NMD is consistent with the well-established EJC model of mammalian NMD. The physical basis for this model was established by the discovery that splicing deposits a multiprotein complex 20–24 nucleotides upstream of exon-exon junctions [1]:

> "This complex protects 8 nucleotides of mRNA from complete RNase digestion at a conserved position 20–24 nucleotides upstream of exon–exon junctions." — Le Hir et al. 2000, *EMBO J*, Abstract

> "In this study, we demonstrate that pre-mRNA splicing stably alters mRNP structure at a conserved position, 20–24 nt upstream of mRNA exon–exon junctions, both in vitro and in vivo." — Le Hir et al. 2000, *EMBO J*, Discussion p.6865

Subsequent work identified RNPS1 as the link between this post-splicing complex and the NMD machinery, through direct interaction with the human Upf complex and a tethered-function experiment demonstrating sufficiency [40]:

> "To ask whether any component of the postsplicing complex constitutes a downstream 'mark' for NMD, we tested the ability of each of the identified proteins to trigger mRNA decay when tethered to the 3′ UTR of β-globin mRNA via fusion to the MS2 coat protein... Of the five postsplicing complex proteins tested, RNPS1 produced a striking down-regulation of β-globin mRNA to 25% of normal... This was comparable to that of the positive control, hUpf3b." — Lykke-Andersen, Shu & Steitz 2001, *Science*, Results p.1837

> "Taken together these results demonstrate that RNPS1, and to a lesser extent Y14, both components of the postsplicing complex, trigger NMD when bound downstream of a translation termination codon." — Lykke-Andersen, Shu & Steitz 2001, *Science*, Results p.1837

> "Our data explain how a dynamic postsplicing complex can function in mRNA quality control... If termination occurs upstream of the last exon-exon junction, interactions between the translation release factors, eRF1 and eRF3, and the downstream postsplicing/hUpf3 complex (via hUpf2 and hUpf1) trigger mRNA decapping followed by rapid decay." — Lykke-Andersen, Shu & Steitz 2001, *Science*, Discussion p.1839

The "50–55 nucleotide rule" — that a stop codon located more than 50–55 nt upstream of an exon-exon junction triggers NMD — was articulated by Nagy and Maquat as a rule drawn from a series of reporter-gene studies [2]:

> "On the basis of studies of transcripts that encode triose phosphate isomerase, the major urinary protein or β-globin, we have defined a rule for termination-codon position: only those termination codons located more than 50–55 nucleotides upstream of the 3′-most exon–exon junction (measured after splicing) mediate a reduction in mRNA abundance." — Nagy & Maquat 1998, *TIBS*, p.198

> "We propose two possible mechanisms by which premature termination of cytoplasmic translation more than 50–55 nucleotides upstream of the 3′-most exon–exon junction causes mRNA decay: a component of the translation-termination complex could interact with a 'mark' positioned at the 3′-most exon–exon junction; alternatively, the cytoplasmic translation-elongation complex might fail to interact with such a mark. We envision association between the mark and mRNA being a consequence of nuclear pre-mRNA splicing." — Nagy & Maquat 1998, *TIBS*, p.198

The mechanistic rationale for the rule — that the terminating ribosome's footprint is insufficient to physically remove a downstream EJC — was later formalized [4, 6]:

> "The 50- to 55-nt rule makes sense considering that a translationally active ribosome poised at a nonsense codon situated more than ∼50–55 nt upstream of an exon–exon junction will not have progressed sufficiently far along the mRNA to remove the EJC deposited ∼20–24 nt upstream of that exon–exon junction. In contrast, it is thought that a ribosome poised at a nonsense codon located either less than ∼50–55 nt upstream of an EJC or downstream from the EJC will have removed the EJC." — Isken & Maquat 2007, *Genes Dev*, pp.1834–1835

> "If assembly of eRF1–eRF3 at a termination codon occurs ≥50–55 nts upstream from an exon-exon junction, the footprint of the terminating ribosome is insufficient to physically remove the EJC, and signal two is engaged." — Popp & Maquat 2013, *Annu Rev Genet*, Step 1: Detection, p.141

The mechanistic pathway from PTC recognition to mRNA degradation proceeds through the SURF complex (SMG1-UPF1-eRF1-eRF3), which assembles at the terminating ribosome and interacts with downstream EJC components to form the DECID complex, committing the transcript to degradation via UPF1 phosphorylation [3, 4]:

> "we describe a novel complex that contains the NMD factors SMG-1 and Upf1, and the translation termination release factors eRF1 and eRF3 (SURF). ... an association between SURF and the EJC is required for SMG-1-mediated Upf1 phosphorylation and NMD." — Kashima et al. 2006, *Genes Dev*, Abstract

> "Based on the essential role of the complex formation of SURF–Upf2–EJC in NMD, we named this complex DECID." — Kashima et al. 2006, *Genes Dev*, Discussion

> "Upf1 binding to the EJC results in the Smg1-mediated phosphorylation of Upf1. Since the Upf factors interact with mRNA degradative activities known to function in NMD, it is reasonable to propose that Upf1 phosphorylation triggers steps that are required for mRNA decay, including the recruitment of degradative activities to mRNA." — Isken & Maquat 2007, *Genes Dev*, p.1842

> "characterizations that include analyses of mRNA decay intermediates indicate that NMD degrades mRNAs primarily by deadenylation-independent decapping via the Dcp1:Dcp2 complex, followed by 5′-to-3′ decay of the transcript body by the Xrn1 exonuclease. ... A minor NMD pathway relies on accelerated deadenylation and, subsequently, 3′-to-5′ decay of the transcript body by the exosome." — Isken & Maquat 2007, *Genes Dev*, p.1843

Our experimental design using SMG1 inhibition directly targets this commitment step.

The EJC-dependent pathway is broadly recognized as the dominant NMD mechanism in mammals [5, 6, 7], though it is not the only one:

> "two different mechanisms have been proposed for target discrimination: 3′ untranslated region (UTR) exon junction complex (EJC)-dependent NMD and 3′ UTR EJC-independent NMD; the latter is also known as long 3′ UTR-mediated NMD or EJC-independent NMD. Generally, a 3′ UTR EJC mediates the most efficient NMD." — Kurosaki, Popp & Maquat 2019, *Nat Rev Mol Cell Biol*, p.408

> "When positioned downstream of a termination codon, the best-characterized role of the EJC is the targeting of mRNAs for and strongly activating NMD... This distinction forms the basis of the '50–55 nucleotide rule': NMD occurs if a PTC is located ≥50–55 nucleotides upstream of an exon–exon junction. In these cases, the leading edge of the terminating ribosome falls short of physically removing the EJC." — Kurosaki, Popp & Maquat 2019, *Nat Rev Mol Cell Biol*, pp.408–409

> "Although the EJC-dependent model for substrate recognition has been rigorously experimentally validated, exceptions exist. NMD is a highly conserved mechanism and the yeast *Saccharomyces cerevisiae*, which has few introns, can recognize aberrant termination on the basis of the distance between the terminating ribosome and the poly(A) tail... Mammalian NMD substrates have also been identified that bear unusually long 3′ UTRs that are devoid of an EJC." — Popp & Maquat 2013, *Annu Rev Genet*, p.145

> "a single 3′UTR EJC that is deposited at an exon–exon junction residing more than ∼50–55 nucleotides downstream of the termination codon... is sufficient to interact with the UPF1 complex, through EJC-bound UPF2 and UPF3 or UPF3X, and trigger NMD. ... [EJC-independent NMD] appears to be less efficient than 3′UTR EJC-promoted NMD." — Kurosaki & Maquat 2016, *J Cell Sci*, p.462

Comparative studies of EJC-dependent and EJC-independent NMD have found the two pathways to be partially redundant, with EJC-enhanced NMD predominating [8]:

> "An exon-junction complex (EJC) located downstream from a TC acts as an NMD-enhancing signal, but is not generally required for NMD. ... our data further indicate the existence of two at least partially redundant decay pathways in NMD of human cells." — Metze et al. 2013, *RNA*, Abstract and Results p.1434

Our model's strong weighting of the `n_downstream_ejc` structural feature (the single most important predictor for PTC+ isoforms) and the dose-response relationship between EJC count and prediction probability quantitatively recapitulate this consensus. Importantly, our model also captures the graded nature of NMD — NMD targets ~10% of the mammalian transcriptome at the isoform level, and the magnitude of NMD-induced downregulation varies by transcript [5, 6, 10]:

> "Although NMD was first found to target one-third of mutated, disease-causing mRNAs, it is now known to also target ~10% of unmutated mammalian mRNAs to facilitate appropriate cellular responses — adaptation, differentiation or death — to environmental changes." — Kurosaki, Popp & Maquat 2019, *Nat Rev Mol Cell Biol*, Abstract

> "In addition to its quality control function, which usually involves mRNA degradation, NMD also controls the abundance of ~10% of the cellular transcriptome. Features that explain why at least some transcripts are recognized by the NMD machinery include an unusually long (>1 kb) 3′ UTR; an uORF with a termination codon that resides ≥50–55 nucleotides upstream of an exon–exon junction; regulated alternative splicing that introduces a PTC; regulated alternative 3′-end formation that results in normal termination codons triggering NMD." — Kurosaki, Popp & Maquat 2019, *Nat Rev Mol Cell Biol*, p.412

> "At least five classes of NMD-inducing features have been described: (a) an upstream open reading frame (uORF) in the 5′ UTR... ; (b) alternative splicing (AS), in which the resulting shift in the translational reading frame generates a stop codon ≥50–55 nts upstream of an exon-exon junction; (c) abnormally long 3′ UTRs; (d) a normal termination codon ≥50–55 nts upstream of an exon-exon junction... ; and (e) UGA codons within certain selenoprotein-encoding mRNAs... The magnitude of NMD-induced downregulation of these unmutated transcripts is generally less than that of authentic PTC-bearing transcripts, leading to the idea that this is a method for fine-tuning rather than eliminating gene expression." — Popp & Maquat 2013, *Annu Rev Genet*, p.147

A clarifying point on the "~10% of the transcriptome" figure is warranted. The reviews above are sometimes read as implying that NMD degrades a substantial fraction of functional, full-length-protein-encoding mRNAs for regulatory purposes. Our data support a narrower interpretation: the ~10% refers to *isoform-level* NMD targets, and these targets are themselves overwhelmingly PTC-bearing isoforms — generated by regulated alternative splicing, alternative 3'-end formation, uORF placement, or intron retention that creates or exposes a PTC. In our isocall dataset, NMD-responsive isoforms are 54–63% PTC-positive by direct annotation; the majority of remaining cases can be attributed to loss of the reference ATG or to downstream-ATG usage that creates a PTC relative to the reference frame (see the PTC-reclassification analysis). The gene-level regulatory consequence of this ~10% is therefore modulation of the *productive-to-NMD isoform ratio*, rather than degradation of mRNAs that would otherwise encode full-length functional protein. The Popp & Maquat enumeration of NMD-inducing features is fully consistent with this interpretation: every class listed (AS-PTC, uORF, long 3'UTR via EJC repositioning, normal-stop-codon-in-penultimate-exon) is a PTC-generating mechanism at the transcript level, even when the host gene continues to produce productive isoforms.

> "Further complexity is added by the discovery that different mRNA substrates require different combinations of NMD trans-effectors. The classical branch of NMD requires the EJC and all of the UPF factors, whereas the fail-safe branch does not require that an EJC reside in the 3′ UTR; however, it does need an EJC situated upstream of the 3′ UTR. Through molecular tethering experiments, a UPF2-independent branch requiring the core EJC components, but not RNPS1, and a pathway requiring UPF2, UPF3X, and RNPS1, but not some of the core EJC components, have been described." — Popp & Maquat 2013, *Annu Rev Genet*, p.145

> "a kinetic competition between efficient translation termination and the assembly of a degradation-triggering NMD complex determines whether an mRNA survives or not." — Karousis, Nasif & Mühlemann 2016, *WIREs RNA*, p.667

---

## 2. 5'UTR-mediated NMD and the "confused ribosome"

Approximately 15% of NMD isoforms in our data lack the canonical PTC + EJC signature, and the model shifts to using start codon context and 5'UTR composition features to predict their NMD status. This subpopulation is consistent with the well-documented class of uORF-containing NMD targets first characterized at scale by Mendell et al. [11], who used microarrays to measure transcript-level responses to NMD inhibition in HeLa cells and found that 4.9% of assayable transcripts were consistently upregulated. uORF-containing mRNAs were one of the major classes identified:

> "We found that 197 transcripts (4.9%) were consistently upregulated and 176 transcripts (4.4%) were consistently downregulated by a factor of at least 1.9 in duplicate experiments." — Mendell et al. 2004, *Nat Genet*, Results p.1075

> "Putative NMD-inducing features included upstream open reading frames (uORFs; 70 transcripts), alternative splicing that introduces nonsense codons or frameshifts (21 transcripts, some of which undergo alternative splicing specifically in HeLa cells) and introns in the 3′ untranslated region (UTR; 9 transcripts)." — Mendell et al. 2004, *Nat Genet*, Results p.1075

> "An intron located at least 50–55 nucleotides downstream of a termination codon is sufficient to initiate mammalian NMD. More than half of the upregulated transcripts (104 of 197) had identifiable features that satisfied this constraint." — Mendell et al. 2004, *Nat Genet*, Results p.1075

Transcriptome-wide UPF1 binding studies confirmed that NMD factors associate with uORF-containing transcripts [12], and that translation is required for the uORF signal:

> "We identified over 200 direct UPF1 binding targets using crosslinking/immunoprecipitation-sequencing (CLIP-seq). ... Translated but not untranslated uORFs are associated with NMD." — Hurt, Robertson & Burge 2013, *Genome Res*

A second line of evidence comes from direct studies of NMD autoregulation: several of the NMD core factors themselves harbour 5'UTR uORFs and are stabilized when NMD is inhibited, demonstrating that uORF-mediated NMD is a bona fide regulatory mechanism rather than an incidental effect [19]:

> "Genes up-regulated in UPF1 and SMG6 knockdown conditions have a higher percentage of such uORF-containing mRNAs. ... UPF2, SMG1, SMG5, and SMG7 mRNAs also contain uORFs, whereas none of the investigated NMD factors harbors an intron >50 nt downstream from the stop codon." — Yepiskoposyan et al. 2011, *RNA*

Multiple reviews have since characterized uORF-mediated NMD as a gene-regulatory mechanism integral to cellular physiology [13, 14, 15]:

> "Stop codons that trigger NMD can also form if the primary coding region of an mRNA is preceded by an upstream open reading frame (uORF). ... NMD also targets non-mutant transcripts, and its regulation of normal gene expression impacts a wide range of physiological processes including cell differentiation, response to stress and development of disease." — Nickless, Bailis & You 2017, *Cell Biosci*

> "The termination codon of a uORF can be recognized as a PTC since it is distant from the 3′UTR signals and the corresponding transcript usually presents downstream EJCs located in the coding sequence of the main ORF." — Barbosa, Peixeiro & Romão 2013, *PLoS Genet*

> "40–50% of human and rodent mRNAs contain at least one uORF. ... uORFs and uAUGs occur less frequently than expected by chance implying that they are under negative selection pressure. ... uORFs reduce protein expression from the downstream mORF by 30–80%." — Somers, Pöyry & Willis 2013, *Int J Biochem Cell Biol*

The mechanistic basis for uORF-mediated NMD involves translation reinitiation. After translating a short uORF, ribosomes can reinitiate at downstream AUGs, but reinitiation efficiency depends on uORF length, inter-cistronic distance, and Kozak context [16]:

> "Figure 1 (lanes 2–4) shows a 3-fold decrease in the efficiency of reinitiation when the upORF was expanded from 13 to 33 codons. ... The critical role of eIF2 was actually predicted, in advance of the GCN4 studies, from the fact that reinitiation became more efficient when the distance between the upORF and the next AUG codon was lengthened." — Kozak 2001, *Nucleic Acids Res*

When reinitiation occurs at a suboptimal position — for example, at an out-of-frame AUG or at a position that places the eventual stop codon upstream of an EJC — the transcript becomes an NMD target. Our model's elevated ATG-branch importance for the PTC- ref ATG retained subgroup, combined with its diffuse attention across multiple ORFs, is consistent with this mechanism: the model searches the 5'UTR landscape for signals of aberrant translation initiation rather than relying on a single canonical start codon.

We describe this pattern informally as the "confused ribosome" hypothesis — that NMD sensitivity in these isoforms arises from uncertainty in ORF selection at the translation initiation stage. While this specific framing is not an established term in the literature, its components are well-supported. The role of Kozak context in modulating start-site selection is foundational, established by Kozak's classic mutagenesis studies [17] and elaborated in the canonical scanning model [18, 39]:

> "From the foregoing mutational analysis, the sequence ACCATGG emerges as the most favorable context for initiation." — Kozak 1986, *Cell*, Discussion p.287

> "Single nucleotide changes in those positions modulate the yield of proinsulin over a 20-fold range. ... The 20-fold variation in proinsulin synthesis among mutants in the B series confirms that sequences flanking the ATG codon modulate translational efficiency." — Kozak 1986, *Cell*, Results pp.285–286

> "The dominant effect of position −3 in eukaryotic ribosome binding sites differs from the prokaryotic Shine-Dalgarno sequence, within which no single position is more important than any other. In striking contrast with the enhanced translation that occurred when A was introduced into position −3, there was little stimulation when A occurred in position −2 or −4." — Kozak 1986, *Cell*, Discussion p.288

> "the optimal context for recognition of the AUG START codon is GCCRCCaugG. Within this motif, the purine (R) in position −3 is the most highly conserved. ... when the first AUG resides in a very weak context, lacking both R in position −3 and G in position +4, some ribosomes initiate at that point but most continue scanning and initiate farther downstream." — Kozak 2002, *Gene*

> "Leaky scanning means that some 40S ribosomal subunits bypass the first AUG codon and initiate instead at the second or, rarely, even the third AUG. ... The most predictable cause of leaky scanning is the absence of a good context around the first AUG codon." — Kozak 1999, *Gene*, §6.1 p.197

> "The fact that many viruses, including HIV, employ this leaky scanning mechanism to produce essential proteins underscores the importance of context (A−3, G+4) in natural situations." — Kozak 1999, *Gene*, §6.1 p.197

The prevalence of uORFs as translational regulators [15, 19] and the regulatory coupling between translation initiation and NMD [14] are all well-documented in the quotes above. The novelty of our contribution lies in the model's data-driven discovery that these features form a coherent alternative prediction strategy, distinct from the canonical PTC pathway, in a measurable subpopulation of NMD targets.

---

## 3. 3'UTR length does not strongly predict NMD

Our analyses provide strong evidence that 3'UTR length does not independently drive NMD susceptibility in our lung cell system. In the isopair analysis (Section 3b of the isoform transition report), three lines of evidence converge: (1) PTC-negative NMD comparators have 3'UTR lengths indistinguishable from Control isoforms (p = 0.37, corrected for PTC-induced 3'UTR inflation by measuring from the reference stop codon position); (2) a prediction model showed that adding 3'UTR length to the full model does not improve holdout AUC (0.896 to 0.897), meaning 3'UTR length carries no predictive information beyond downstream EJC count; and (3) 3'UTR splicing is depleted in PTC-negative NMD comparators relative to Control (OR = 0.71, p = 2.63e-8), while enriched in PTC-positive comparators (OR = 1.78), consistent with 3'UTR splicing operating as a PTC-creating mechanism (via EJC repositioning) rather than an independent NMD trigger. Among same-stop PTC+ pairs — where the comparator shares the reference stop codon but has gained 3'UTR splice junctions — 86% have more downstream EJCs in the comparator, confirming the EJC repositioning mechanism.

These findings engage the "faux 3'UTR" model, originally proposed by Amrani et al. [20] in yeast, which posits that when a stop codon is far from the poly(A) tail, the resulting 3'UTR is perceived as aberrant and triggers NMD. The yeast experiments directly tested this by tethering PABP downstream of a PTC:

> "Recent experiments challenge this notion and suggest a model that posits that mRNA decay is activated by the intrinsically aberrant nature of premature termination. Here we use a primer extension inhibition (toeprinting) assay to delineate ribosome positioning and find that premature translation termination in yeast extracts is indeed aberrant. Ribosomes encountering premature UAA or UGA codons in the *CAN1* mRNA fail to release and, after release, migrate to initiation AUGs. This anomaly depends on prior nonsense codon recognition and is eliminated in extracts lacking the principal NMD factor, Upf1p, or by flanking the nonsense codon with a normal 3′-untranslated region (UTR)." — Amrani et al. 2004, *Nature*, Abstract/opening p.112

> "Pab1p tethered 37–73 nt 3′ to premature UAA, UGA or UAG codons promoted 5–11-fold increases in mRNA stability and abundance. As controls, MS2 dimer or MS2–Sxl (sex lethal protein) tethered at the same position, or MS2–Pab1p tethered 164 nt downstream of the stop codon and 3′ to the DSE, had no effect on mRNA stability or abundance." — Amrani et al. 2004, *Nature*, Results p.115

> "Our results indicate that not all termination events are equivalent and that, at least in yeast, NMD is triggered by a ribosome's failure to terminate adjacent to a properly configured 3′-UTR. As suggested by the *faux* UTR model, proper termination of translation and normal rates of mRNA decay are likely to require interactions between a terminating ribosome and a specific RNP structure or set of factors (including Pab1p and Sup35p) localized 3′ to the stop codon." — Amrani et al. 2004, *Nature*, Discussion p.117

In yeast, which lacks EJCs, this distance-dependent mechanism is the primary NMD trigger, and direct evidence shows that wild-type yeast transcripts with long 3'UTRs are preferentially degraded by NMD [35]:

> "Long 3′-UTRs target wild-type mRNAs for decay by the NMD pathway. ... Wild-type yeast mRNAs with exceptionally long 3′-UTRs are strongly enriched for degradation by the NMD pathway. ... Replacement of a long 3′-UTR with a shorter 3′-UTR is sufficient to prevent degradation of the wild-type PGA1 mRNA by NMD." — Kebaara & Atkin 2009, *Nucleic Acids Res*

In mammals, the picture is more nuanced. Several studies have demonstrated that long 3'UTRs can trigger NMD in mammalian cells [21], with the Ig-μ system providing a direct test:

> "If the same mechanism exists in mammals, it should be possible to convert our PTC-free mini-μ WT mRNA into an NMD substrate by extending the 3′ UTR. Indeed, a construct termed mini-μ long 3′ UTR, which has inserted into the 3′ UTR a 1.3-kilobase stuffer sequence, gave a five-fold lower mRNA level compared to the parental mini-μ construct. Knockdown of human Upf1 using RNAi techniques resulted in an increase of mini-μ long 3′ UTR mRNA by four-fold... Thus, increasing the distance between the termination codon and the poly(A) tail in mini-μ leads to a Upf1-dependent mRNA downregulation." — Bühler et al. 2006, *Nat Struct Mol Biol*, Results pp.463–464

> "the results presented here are incompatible with the current model for mammalian NMD, according to which recognition of a PTC strictly relies on an interaction between the terminating ribosome and an EJC located downstream of it. ... our results suggest an NMD model where PTC recognition is EJC independent, evolutionarily conserved... and where the presence of an EJC downstream of the termination codon functions as an enhancer of NMD in mammals." — Bühler et al. 2006, *Nat Struct Mol Biol*, Discussion p.464

UPF1 itself has been shown to bind 3'UTRs in a length-dependent manner [22]:

> "Upf1 copurification with tagged mRNAs increased with 3′UTR length. The relationship between Upf1 copurification and 3′UTR length was strikingly linear, consistent with sequence-nonspecific recognition of long 3′UTRs by Upf1." — Hogg & Goff 2010, *Cell*, Results p.382

> "these findings indicate that the extent of Upf1 association with a transcript is diagnostic of its NMD susceptibility, consistent with previous experiments in yeast, *C. elegans*, and human cells." — Hogg & Goff 2010, *Cell*, Results p.383

> "Based on our findings that rare readthrough events permit Upf1 3′UTR length-dependent accumulation in mRNPs but inhibit decay, we propose a two-step model in which Upf1 senses 3′UTR length to potentiate decay. In this model, length-dependent equilibrium binding of Upf1 marks 3′UTRs of potential decay substrates and increases the probability of Upf1 binding to release factors." — Hogg & Goff 2010, *Cell*, Discussion p.386

The competition model proposes that NMD is determined by the balance between UPF1/EJC complexes (pro-NMD) and PABP proximity to the terminating ribosome (anti-NMD), with longer 3'UTRs tipping the balance toward degradation [23, 24]:

> "the proximity of PABPC1 provides an important signal for defining a translation termination event as 'correct' and prevents degradation of the mRNA by NMD. ... the physical distance between the TC and the poly(A) tail might be a crucial determinant to identify a TC as premature. ... the two antagonizing signals (PABPC1 proximity and Upf1–3 recruitment, respectively) determine if a translation termination event is defined as premature or correct. ... the NMD-inhibiting signal (poly(A) tail proximity) efficiently competes with the NMD-promoting signal (the downstream EJC)." — Eberle et al. 2008, *PLoS Biol*

> "PTC recognition is determined by a competition between 3′ UTR–associated factors, which stimulate (including the EJC) or antagonize (including cytoplasmic PABP) the recruitment of the Upf complex to the terminating ribosome. ... a translation termination event proximal to cytoplasmic PABP, or other unknown NMD-antagonizing factors, precludes the interaction of hUpf1 with eRF3. ... if hUpf1 associates with eRF3, NMD ensues. This occurs when cytoplasmic PABP, or other inhibitory factors, are spatially distant from the termination event." — Singh, Rebbapragada & Lykke-Andersen 2008, *PLoS Biol*

However, several lines of evidence suggest this pathway plays a limited role in mammalian steady-state NMD. Many long 3'UTRs contain protective elements — including structured RNA elements and PTBP1 binding sites — that actively prevent NMD [25, 26]:

> "many [long 3'UTR transcripts] are indeed resistant to NMD. ... [we] identified a cis element located within the first 200 nt that inhibits NMD. ... a subset of long 3′ UTR mRNAs evades NMD by a different mechanism. ... NMD-evading cis elements could also provide an additional level of post-transcriptional gene regulation." — Toma et al. 2015, *RNA*

> "When bound near a stop codon, PTBP1 blocks the NMD protein UPF1 from binding 3'UTRs. ... PTBP1 functions to exclude UPF1 from 3'UTRs, disrupting its ability to accurately discriminate 3'UTR length and induce decay. ... PTBP1 assembled on multiple TC-proximal CU-rich sequences could thus result in a zone of UPF1 exclusion." — Ge, Quek, Beemon & Hogg 2016, *eLife*

A growing consensus holds that the faux 3'UTR model is more relevant in yeast than in mammals, where the EJC pathway dominates — with the WIREs RNA review by Kishor, Fritz & Hogg summarizing the current synthesis [27]:

> "To explain these results, a 'faux 3′UTR' model posits that increased distance between the terminating ribosome and the poly-A tail disfavors the interaction between PABP and eRF3, impairing termination efficiency. In combination with enhanced recognition of long 3′UTRs by UPF1, this model may help explain why the probability of decay increases with 3′UTR length." — Kishor, Fritz & Hogg 2019, *WIREs RNA*, §5.1.3 p.9

> "The hypothesis that NMD recognizes termination events due to the absence of a functional PABP-eRF3 interaction has been challenged by efforts in yeast to directly test the model: elimination of PABP-eRF3 binding, transient depletion of PABP, and use of reporter mRNAs lacking poly-A tails all failed to generate the predicted NMD defects." — Kishor, Fritz & Hogg 2019, *WIREs RNA*, §5.1.3 p.9

> "In complex vertebrate transcriptomes, the EJC is a molecular marker that enhances the discrimination of NMD substrates. Transcriptomes in eukaryotes that use the EJC for NMD have coevolved with the pathway such that it is relatively rare for stop codons to occur more than 50 nt upstream of the final exon-exon junction. Given that EJCs are deposited ~20 nucleotides upstream of an exon-exon boundary, EJCs remaining more than ~30 nucleotides downstream of the TC can thus be used as a marker of premature termination, providing a rationale for degradation by the NMD pathway." — Kishor, Fritz & Hogg 2019, *WIREs RNA*, §5.2 p.9

> "Mechanistic flexibility may be another strategy evolved by the NMD pathway to function effectively in complex transcriptomes, since the streamlined yeast pathway shows equal dependence on UPF1, UPF2, and UPF3 for degradation of NMD target mRNAs." — Kishor, Fritz & Hogg 2019, *WIREs RNA*, §6 p.13

Our experimental data therefore align well with the emerging view: in the context of endogenous human isoform diversity, 3'UTR length is not an independent NMD trigger. When 3'UTR splicing does contribute to NMD, it does so by repositioning EJCs to create PTCs — reinforcing rather than replacing the canonical EJC-dependent pathway.

---

## 4. Stochastic ORF selection and NMD susceptibility

We hypothesize that ribosomal ORF selection may be partially stochastic, such that NMD susceptibility of a transcript is a function of start codon strength (Kozak context) combined with the "NMD potential" of the downstream reading frame. This hypothesis is supported by our model's multi-ORF architecture, in which the learned attention mechanism distributes weight across up to 5 candidate ORFs per transcript — particularly for non-canonical NMD isoforms where attention entropy is high.

The mechanistic foundation for stochastic ORF selection is the leaky scanning model: the 43S preinitiation complex can bypass AUGs in suboptimal Kozak context, allowing ribosomes to initiate at downstream start codons [17, 18, 28]:

> "most mRNAs are translated by a scanning mechanism wherein the small (40S) ribosomal subunit is preloaded with Met-tRNAi by the GTP-bound form of eukaryotic initiation factor 2 (eIF2)... the resulting 43S preinitiation complex (PIC) then attaches to the mRNA... Attachment of the 43S complex is confined to the free 5′ end of the mRNA... and the 5′ untranslated region (5′UTR) is scanned base by base for complementarity to the anticodon (AC) of Met-tRNAi as successive triplets enter the P site of the 40S subunit. Thus, the first AUG encountered is favored as the start codon." — Hinnebusch 2014, *Annu Rev Biochem*, Overview p.780

> "particular sequences immediately surrounding the AUG, especially those including a purine at position −3, enhance AUG selection by the scanning PIC; and a 5′-proximal AUG that deviates sufficiently from the optimum context, which in mammals is 5′-(A/G)NNAUGG-3′, can be bypassed in an event termed **leaky scanning**. Shortening the 5′UTR beyond ~20 nt also reduces the efficiency of initiation, a finding that can be exploited to produce an N-terminally extended polypeptide (by inefficient initiation at the 5′-proximal AUG) in addition to the shorter, major isoform (by efficient initiation at the downstream AUG)." — Hinnebusch 2014, *Annu Rev Biochem*, Overview p.780

> "A critical aspect of the scanning process is the ability of the 43S PIC to bypass AUGs in poor surrounding sequence context, as well as near-cognate triplets (those with single-base mismatches from AUG) in the 5′UTR, so that the initiation complex can be assembled at the correct AUG start codon on the mRNA." — Hinnebusch 2014, *Annu Rev Biochem*, AUG Recognition by the Scanning PIC p.793

> "The nature of scanning, its 5′ to 3′ directionality, dictates that the initiation codon is frequently the AUG triplet closest to the 5′ end, encountered first by the scanning PIC. The first AUG can be skipped when it is flanked by an unfavorable sequence—a process termed 'leaky scanning'—to use a downstream AUG. A favorable sequence context in mammals is the 'Kozak consensus', 5′ (A/G)CCAUGG 3′." — Hinnebusch, Ivanov & Sonenberg 2016, *Science*, pp.1–2

> "If the uAUG is followed by a stop codon in the same ORF, then translation of the upstream ORF (uORF) will attenuate translation of the downstream main ORF, because reinitiation is generally inefficient." — Hinnebusch, Ivanov & Sonenberg 2016, *Science*, p.2

> "Termination at an uORF stop codon can elicit the same mRNA destabilization evoked by the nonsense-mediated decay (NMD) pathway at premature termination codons in ORFs, magnifying the inhibitory effects of uORFs." — Hinnebusch, Ivanov & Sonenberg 2016, *Science*, p.4

> "uORFs whose AUG codons better conform to the Kozak consensus are more inhibitory... upstream start codons tend to be near-cognates or AUGs in poor context, which should favor leaky-scanning." — Hinnebusch, Ivanov & Sonenberg 2016, *Science*, Translational control by uORFs p.4

Ribosome profiling studies have provided direct evidence that single transcripts engage ribosomes at multiple positions. The foundational work by Ingolia et al. [29, 30] revealed widespread translation of uORFs and alternative ORFs:

> "The position of a translating ribosome can be precisely determined using the fact that a ribosome protects a discrete footprint [~30 nucleotides] on its mRNA template from nuclease digestion... Here, we present a ribosome profiling strategy based on deep sequencing of ribosome-protected fragments that provides comprehensive, high precision measurements of in vivo translation with sub-codon precision." — Ingolia et al. 2009, *Science*, Introduction

> "Comparing the rate of translation to mRNA abundance from the same samples revealed a roughly 100-fold range of translation efficiency... between different yeast genes... Thus, differences in translational efficiency, which are invisible to mRNA abundance measurements, contribute substantially to the dynamic range of gene expression." — Ingolia et al. 2009, *Science*, Results p.2

> "More broadly, among all annotated 5′ UTRs we found evidence for the translation of 153 uORFs, fewer than 30 of which had been experimentally evaluated previously. ... Overall, we found 143 non-AUG uORFs with evidence of translation, which account for 20% of 5′UTR ribosome footprints. Thus, there is pervasive initiation at specific, favorable non-AUG sites." — Ingolia et al. 2009, *Science*, Results p.3

> "We find thousands of pauses in the body of genes (1500 pauses in 1100 genes in a set of 4994 well-expressed genes) and at termination codons (420 pauses)." — Ingolia, Lareau & Weissman 2011, *Cell*, Results p.791

> "The majority of unannotated near-cognate initiation sites we detected drive the translation of upstream open reading frames (uORFs)... These uORF initiation sites are accompanied by elongating ribosome footprints in the untreated sample that are depleted during harringtonine treatment, indicating that they are involved in active translation... Our observations suggest that near-cognate uORFs are quite common." — Ingolia, Lareau & Weissman 2011, *Cell*, Results p.796

> "A number of features of mammalian proteomes emerge from our studies, including the ubiquitous use of alternate initiation sites that drive the production of extended or truncated isoforms of known proteins as well as the translation of sprcRNAs, whose protein-coding potential was not initially apparent. We also observe widespread translation upstream of mammalian protein-coding genes, similar to but more extensive than upstream translation that we observed in yeast." — Ingolia, Lareau & Weissman 2011, *Cell*, Discussion p.799

Specialized initiation-site profiling confirmed that most mRNAs have multiple active start codons with usage frequencies that correlate with — but are not fully determined by — Kozak context strength [31, 32]:

> "Remarkably, nearly half the transcripts (49.6%) contained multiple TIS sites, suggesting that alternative translation prevails even under physiological conditions. ... GTI-seq revealed that 54% of transcripts bear one or more TIS positions upstream of the annotated start codon. ... The AUG codon with a strong Kozak sequence context showed higher initiation efficiency (or lower leakiness) than a codon with a weak or no consensus sequence. ... ribosome leaky scanning tends to occur when the context for an aTIS is suboptimal." — Lee et al. 2012, *PNAS*

> "many starvation-responsive genes contain multiple TISs (1,286 in HEK293 and 1,343 in MEF), suggesting a regulatory role for alternative TISs in translational control. ... Among transcripts with increased aTIS initiation upon starvation, the Kozak consensus motif is prominent." — Gao et al. 2015, *Nat Methods*

Cross-species analyses have shown that uORF-mediated translational repression is conserved and quantitatively related to Kozak context [33, 34]:

> "the regulatory potential of uORFs on individual genes is conserved across species. ... uORF initiation contexts display signatures of selection. ... High- and medium-confidence translated uORFs possess significantly better initiation contexts than background." — Johnstone, Bazzini & Giraldez 2016, *EMBO J*

> "the repressiveness and sequence features of uORFs are conserved in vertebrates. ... the ratio of translation over 5′ leaders and CDSes is conserved between human and mouse, and correlates with the number of uORFs. ... more favourable initiation context sequences and less-stable secondary structures correlate with increased uORF TE." — Chew, Pauli & Schier 2016, *Nat Commun*

Synthetic-biology studies further demonstrate that tuning the translation-initiation context of a uORF quantitatively controls main-ORF output [38]:

> "By varying the base sequence preceding the uORF, we sought to vary the translation initiation rate of the uORF and subsequently control the degree of this suppression. ... The uORF acts to shunt ribosomes so that fewer ribosomes reach the downstream ORF. ... By tuning the TIS of both a single uORF and a downstream primary ORF, we were able to achieve a nearly continuous range of translation from 0.05 to 0.6 relative units." — Ferreira, Overton & Wang 2013, *PNAS*

The direct connection between ribosome occupancy at uORF termination codons and NMD susceptibility was first demonstrated in the yeast CPA1 system, where a peptide-dependent stall at the uORF stop codon is required to trigger NMD [36]:

> "Our data indicate that ribosome stalling at the uORF termination codon induces NMD. Consistent with this, the *CPA1* mRNA... decays more rapidly in NMD⁺ cells grown in medium containing Arg than in medium lacking Arg. ... Arg addition to media rapidly destabilized the *CPA1* transcript in wild-type but not *upf1Δ* cells." — Gaba, Jacobson & Sachs 2005, *Mol Cell*, Summary p.449 and Discussion p.456

> "The native context wild-type CPA1 uORF efficiently stalls ribosomes and triggers NMD, while the D13N uORF, which stalls ribosomes poorly or not at all, is an inefficient trigger. Improving the uORF initiation context renders both the wild-type and D13N uORFs efficient triggers of NMD. The NMD-inducing effects of the context-improved D13N uORF were not due to ribosome stalling because the D13N uORF did not mediate ribosome-stalling activity, regardless of its initiation context. These data indicate that the levels of ribosomes at an early stop codon can be regulated to control NMD." — Gaba, Jacobson & Sachs 2005, *Mol Cell*, Discussion p.457

> "These experiments showed that CPA1 mRNA half-lives were ~7 and 3 min, respectively, for NMD⁺ cells grown in media that lacked or contained Arg, while, in upf1Δ cells, the corresponding half-lives were ~14 and 17 min, respectively." — Gaba, Jacobson & Sachs 2005, *Mol Cell*, Results p.452

Our Kozak PWM analysis supports this framework: PTC+ isoforms have the strongest Kozak context at their reference CDS start codons (mean PWM = 0.79), while PTC- ref ATG retained isoforms — despite having the same reference CDS ORF — show weaker Kozak context (mean PWM = 0.30), potentially allowing more ribosomal read-through to alternative ORFs with different NMD fates.

The specific formulation that NMD susceptibility should be modeled as a probability-weighted function over multiple candidate ORFs — where the probability is governed by Kozak context strength and the NMD potential by the ORF's termination context — appears to be a novel analytical framework. While the individual components are well-established, we are not aware of prior work that has explicitly modeled this interaction in a multi-ORF deep learning architecture or demonstrated its predictive utility at transcriptome scale.

---

## 5. Are there genuinely non-PTC NMD mechanisms in mammals?

A recurring source of confusion in the NMD literature is the conflation of two different definitions of "PTC":

- **Pathological PTC** — a stop codon introduced by a disease-causing mutation upstream of the natural stop codon.
- **Operational PTC** — any stop codon that sits ≥50–55 nt upstream of a downstream exon-exon junction (the Nagy–Maquat rule [2]).

Under the operational definition, almost all mammalian NMD is PTC-based. Most of the mechanisms commonly described as "non-PTC" or "regulatory" NMD turn out to be operational-PTC NMD under a different name:

- **uORF-mediated NMD** — the uORF termination codon is itself an operational PTC relative to the main ORF; the EJCs that trigger decay sit downstream of the uORF stop, in the coding sequence of the main ORF [14].
- **AS-NMD (regulated alternative splicing introducing a PTC)** — by definition creates a PTC in the alternative isoform [6].
- **3'UTR introns** (Mendell's 9 transcripts, ~5% of NMD-regulated transcripts [11]) — the stop codon is the *normal* one, but the downstream intron deposits an EJC >50 nt past the stop, making the normal stop an operational PTC.
- **Selenoprotein UGA codons** — UGA is a conditional PTC that acts as a stop when Se is depleted.

This leaves a narrow window for *genuinely* PTC-independent NMD — mRNAs where the stop codon has no downstream EJC yet the transcript is still an NMD substrate. The literature offers three candidate mechanisms:

### 5.1 Long 3'UTR-mediated, EJC-independent NMD

This is the strongest candidate. Bühler et al. [21] took a PTC-free mini-μ construct with no downstream intron and made it an NMD substrate purely by extending the 3'UTR:

> "a construct termed mini-μ long 3′ UTR, which has inserted into the 3′ UTR a 1.3-kilobase stuffer sequence, gave a five-fold lower mRNA level compared to the parental mini-μ construct. Knockdown of human Upf1 using RNAi techniques resulted in an increase of mini-μ long 3′ UTR mRNA by four-fold... Thus, increasing the distance between the termination codon and the poly(A) tail in mini-μ leads to a Upf1-dependent mRNA downregulation." — Bühler et al. 2006, *Nat Struct Mol Biol*, Results pp.463–464

Hogg & Goff [22] provided the mechanistic backbone: UPF1 binds 3'UTRs in a length-dependent, sequence-nonspecific manner, and UPF1 coverage correlates with NMD susceptibility. However, this mechanism has three serious limitations in endogenous mammalian transcriptomes:

1. **Protective cis-elements neutralize long 3'UTRs.** Toma et al. [25] and Ge et al. [26] show that many endogenous long 3'UTRs are NMD-resistant because they contain structured RNA elements and PTBP1 binding sites that exclude UPF1.

2. **Direct tests of the faux-UTR model have failed in yeast.** Kishor, Fritz & Hogg [27] summarize the current state:
> "The hypothesis that NMD recognizes termination events due to the absence of a functional PABP-eRF3 interaction has been challenged by efforts in yeast to directly test the model: elimination of PABP-eRF3 binding, transient depletion of PABP, and use of reporter mRNAs lacking poly-A tails all failed to generate the predicted NMD defects." — Kishor, Fritz & Hogg 2019, *WIREs RNA*, §5.1.3 p.9

3. **Our isopair data argue directly against this mechanism in endogenous human isoforms.** 3'UTR length, measured from the reference stop codon position (so as not to conflate with PTC-induced inflation), is indistinguishable between PTC-negative NMD comparators and Controls (p = 0.37). In a prediction model, adding 3'UTR length provides essentially no AUC improvement over downstream EJC count alone (0.896 → 0.897). If the long-3'UTR mechanism operated at endogenous scale, both tests should have been positive.

### 5.2 UPF1-CLIP "non-canonical" targets (Hurt 2013)

In their genome-wide UPF1 CLIP-seq study, Hurt, Robertson & Burge [12] describe a subset of UPF1-repressed transcripts that do not fit either the PTC+downstream-EJC or the long-3'UTR model:

> "A repression pathway that involves 3′ UTR binding by UPF1 and translation but is independent of canonical targeting features involving 3′ UTR length and stop codon placement." — Hurt, Robertson & Burge 2013, *Genome Res*

The paper identifies this pathway by exclusion (UPF1-bound + translation-dependent + lacks canonical features) but does not mechanistically characterize it. To our knowledge, no subsequent study has converged on a concrete molecular model for this population.

### 5.3 Yeast long-3'UTR NMD

In *S. cerevisiae*, which lacks EJCs, long-3'UTR NMD is the dominant mechanism [35]. This is a genuinely non-PTC pathway in yeast — but it does not transfer cleanly to mammals, where the EJC pathway predominates and where direct tests of the analogous mammalian mechanism have been equivocal (§5.1).

### Synthesis

The most defensible reading of the current literature is that **mammalian NMD is overwhelmingly PTC-based under the operational definition**. The reviews that report NMD regulating ~10% of the "normal" transcriptome are correctly read as: NMD targets ~10% of expressed isoforms, and those isoforms are themselves predominantly operational-PTC-bearing, generated by regulated splicing, alternative 3'-end formation, uORF placement, or intron retention. The gene-level regulatory consequence arises from shifting the productive-to-NMD isoform ratio, not from NMD of functional full-length-protein-encoding mRNAs.

The one candidate mechanism that would genuinely require a non-PTC model — long-3'UTR / UPF1-density-driven decay — appears to be narrow in scope, reporter-dependent, and often neutralized by protective cis-elements in endogenous transcripts. Our isopair data provide two orthogonal tests that this mechanism is weak-to-absent in endogenous human isoforms, supporting the narrower interpretation.

---

## Citation corrections

During quote extraction, four entries in the original bibliography were found to be incorrect. All four have been corrected in the reference list below; the corrections are itemized here for transparency.

1. **Ref 9 (originally "Bhatt DM et al. 2012, *Mol Cell* 46(5):585-595, 'Transcript dynamics of proinflammatory genes reveal that p65 mediates a revised NMD pathway'")** — **no such paper exists**. The closest matching real paper is Bhatt et al. 2012 *Cell* 150(2):279-290 ("Transcript Dynamics of Proinflammatory Genes Revealed by Sequence Analysis of Subcellular RNA Fractions"), which does not discuss NMD. The original Ref 9 has been removed from the bibliography; the graded-NMD claim now relies on Karousis, Nasif & Mühlemann 2016 (Ref 10), Kurosaki et al. 2019 (Ref 5), and Popp & Maquat 2013 (Ref 6).

2. **Ref 10 (journal error)** — the paper "Nonsense-mediated mRNA decay: novel mechanistic insights and biological impact" by Karousis, Nasif & Mühlemann was published in **WIREs RNA 7(5):661-682 (2016)**, not *Biol Chem* 397:1093-1126. The PMC ID is PMC6680220. The citation has been corrected.

3. **Ref 27 (wrong first author/year)** — the paper titled "Nonsense-mediated mRNA decay: the challenge of telling right from wrong in a complex transcriptome" was authored by **Kishor A, Fritz SE, Hogg JR** and published in **WIREs RNA 10(6):e1548 (2019)**, not by Boehm/Haberman/Hentze/Kulozik in 2014. The citation has been corrected.

4. **Ref 35 (Kebaara & Atkin 2009)** — this paper is about *3'UTR length driving NMD in yeast* and was incorrectly cited in the original Section 4 for the claim "altering translation initiation efficiency affects NMD susceptibility." The Kebaara citation has been moved to Section 3 where it correctly supports the yeast faux-3'UTR pathway. In Section 4, the translation-initiation-efficiency claim now relies on Gaba, Jacobson & Sachs 2005 (Ref 36) alone.

---

## Papers folder

PDFs for all 39 valid references are in `./papers/`. Files follow the pattern `ref{N}_FirstAuthorYear.pdf` (with minor variations where the user placed papers with slightly different naming). Every quoted passage in this review has been verified against the downloaded PDF.

---

## References

1. Le Hir H, Izaurralde E, Maquat LE, Moore MJ. The spliceosome deposits multiple proteins 20-24 nucleotides upstream of mRNA exon-exon junctions. *EMBO J.* 2000;19(24):6860-6869. [PMC305905](https://pmc.ncbi.nlm.nih.gov/articles/PMC305905/)

2. Nagy E, Maquat LE. A rule for termination-codon position within intron-containing genes: when nonsense affects RNA abundance. *Trends Biochem Sci.* 1998;23(6):198-199. [PubMed 9644970](https://pubmed.ncbi.nlm.nih.gov/9644970/)

3. Kashima I, Yamashita A, Izumi N, et al. Binding of a novel SMG-1-Upf1-eRF1-eRF3 complex (SURF) to the exon junction complex triggers Upf1 phosphorylation and nonsense-mediated mRNA decay. *Genes Dev.* 2006;20(3):355-367. [PMC1361706](https://pmc.ncbi.nlm.nih.gov/articles/PMC1361706/)

4. Isken O, Maquat LE. Quality control of eukaryotic mRNA: safeguarding cells from abnormal mRNA function. *Genes Dev.* 2007;21(15):1833-1856. [PubMed 17671086](https://pubmed.ncbi.nlm.nih.gov/17671086/)

5. Kurosaki T, Popp MW, Maquat LE. Quality and quantity control of gene expression by nonsense-mediated mRNA decay. *Nat Rev Mol Cell Biol.* 2019;20(7):406-420. [PMC6855384](https://pmc.ncbi.nlm.nih.gov/articles/PMC6855384/)

6. Popp MW, Maquat LE. Organizing principles of mammalian nonsense-mediated mRNA decay. *Annu Rev Genet.* 2013;47:139-165. [PMC4148824](https://pmc.ncbi.nlm.nih.gov/articles/PMC4148824/)

7. Kurosaki T, Maquat LE. Nonsense-mediated mRNA decay in humans at a glance. *J Cell Sci.* 2016;129(3):461-467. [PMC4760306](https://pmc.ncbi.nlm.nih.gov/articles/PMC4760306/)

8. Metze S, Herzog VA, Ruepp MD, Mühlemann O. Comparison of EJC-enhanced and EJC-independent NMD in human cells reveals two partially redundant degradation pathways. *RNA.* 2013;19(10):1432-1448. [PMC3854533](https://pmc.ncbi.nlm.nih.gov/articles/PMC3854533/)

9. *[removed — original citation was fabricated; see Citation corrections]*

10. Karousis ED, Nasif S, Mühlemann O. Nonsense-mediated mRNA decay: novel mechanistic insights and biological impact. *Wiley Interdiscip Rev RNA.* 2016;7(5):661-682. [PMC6680220](https://pmc.ncbi.nlm.nih.gov/articles/PMC6680220/) *(journal corrected from original; see Citation corrections)*

11. Mendell JT, Sharifi NA, Meyers JL, Martinez-Murillo F, Dietz HC. Nonsense surveillance regulates expression of diverse classes of mammalian transcripts and mutes genomic noise. *Nat Genet.* 2004;36(10):1073-1078. [PubMed 15448691](https://pubmed.ncbi.nlm.nih.gov/15448691/)

12. Hurt JA, Robertson AD, Burge CB. Global analyses of UPF1 binding and function reveal expanded scope of nonsense-mediated mRNA decay. *Genome Res.* 2013;23(10):1636-1650. [PMC3787261](https://pmc.ncbi.nlm.nih.gov/articles/PMC3787261/)

13. Nickless A, Bailis JM, You Z. Control of gene expression through the nonsense-mediated RNA decay pathway. *Cell Biosci.* 2017;7:26. [PMC5437625](https://pmc.ncbi.nlm.nih.gov/articles/PMC5437625/)

14. Barbosa C, Peixeiro I, Romão L. Gene expression regulation by upstream open reading frames and human disease. *PLoS Genet.* 2013;9(8):e1003529. [PMC3738444](https://pmc.ncbi.nlm.nih.gov/articles/PMC3738444/)

15. Somers J, Pöyry T, Willis AE. A perspective on mammalian upstream open reading frame function. *Int J Biochem Cell Biol.* 2013;45(8):1690-1700. [PMC7172355](https://pmc.ncbi.nlm.nih.gov/articles/PMC7172355/)

16. Kozak M. Constraints on reinitiation of translation in mammals. *Nucleic Acids Res.* 2001;29(24):5226-5232. [PMC97554](https://pmc.ncbi.nlm.nih.gov/articles/PMC97554/)

17. Kozak M. Point mutations define a sequence flanking the AUG initiator codon that modulates translation by eukaryotic ribosomes. *Cell.* 1986;44(2):283-292. [PubMed 3943125](https://pubmed.ncbi.nlm.nih.gov/3943125/)

18. Kozak M. Pushing the limits of the scanning mechanism for initiation of translation. *Gene.* 2002;299(1-2):1-34. [PMC7126118](https://pmc.ncbi.nlm.nih.gov/articles/PMC7126118/)

19. Yepiskoposyan H, Aeschimann F, Nilsson D, Okoniewski M, Mühlemann O. Autoregulation of the nonsense-mediated mRNA decay pathway in human cells. *RNA.* 2011;17(12):2108-2118. [PMC3222124](https://pmc.ncbi.nlm.nih.gov/articles/PMC3222124/)

20. Amrani N, Ganesan R, Kerber S, Ghosh S, Jacobson A. A faux 3'-UTR promotes aberrant termination and triggers nonsense-mediated mRNA decay. *Nature.* 2004;432(7013):112-118. [PubMed 15525991](https://pubmed.ncbi.nlm.nih.gov/15525991/)

21. Bühler M, Steiner S, Mohn F, Paillusson A, Mühlemann O. EJC-independent degradation of nonsense immunoglobulin-mu mRNA depends on 3'UTR length. *Nat Struct Mol Biol.* 2006;13(5):462-464. [PubMed 16622410](https://pubmed.ncbi.nlm.nih.gov/16622410/)

22. Hogg JR, Goff SP. Upf1 senses 3'UTR length to potentiate mRNA decay. *Cell.* 2010;143(3):379-389. [PMC2981159](https://pmc.ncbi.nlm.nih.gov/articles/PMC2981159/)

23. Eberle AB, Stalder L, Mathys H, Orozco RZ, Mühlemann O. Posttranscriptional gene regulation by spatial rearrangement of the 3' untranslated region. *PLoS Biol.* 2008;6(4):e92. [PMC2689704](https://pmc.ncbi.nlm.nih.gov/articles/PMC2689704/)

24. Singh G, Rebbapragada I, Lykke-Andersen J. A competition between stimulators and antagonists of Upf complex recruitment governs human nonsense-mediated mRNA decay. *PLoS Biol.* 2008;6(4):e111. [PMC2689706](https://pmc.ncbi.nlm.nih.gov/articles/PMC2689706/)

25. Toma KG, Rebbapragada I, Durand S, Lykke-Andersen J. Identification of elements in human long 3'UTRs that inhibit nonsense-mediated decay. *RNA.* 2015;21(5):887-897. [PMC4408796](https://pmc.ncbi.nlm.nih.gov/articles/PMC4408796/)

26. Ge Z, Quek BL, Beemon KL, Hogg JR. Polypyrimidine tract binding protein 1 protects mRNAs from recognition by the nonsense-mediated mRNA decay pathway. *eLife.* 2016;5:e11155. [PMC4764554](https://pmc.ncbi.nlm.nih.gov/articles/PMC4764554/)

27. Kishor A, Fritz SE, Hogg JR. Nonsense-mediated mRNA decay: the challenge of telling right from wrong in a complex transcriptome. *Wiley Interdiscip Rev RNA.* 2019;10(6):e1548. [PMC6788943](https://pmc.ncbi.nlm.nih.gov/articles/PMC6788943/) *(corrected from original; see Citation corrections)*

28. Hinnebusch AG. The scanning mechanism of eukaryotic translation initiation. *Annu Rev Biochem.* 2014;83:779-812. [PubMed 24499181](https://pubmed.ncbi.nlm.nih.gov/24499181/)

29. Ingolia NT, Ghaemmaghami S, Newman JRS, Weissman JS. Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling. *Science.* 2009;324(5924):218-223. [PMC2746483](https://pmc.ncbi.nlm.nih.gov/articles/PMC2746483/)

30. Ingolia NT, Lareau LF, Weissman JS. Ribosome profiling of mouse embryonic stem cells reveals the complexity and dynamics of mammalian proteomes. *Cell.* 2011;147(4):789-802. [PMC3225288](https://pmc.ncbi.nlm.nih.gov/articles/PMC3225288/)

31. Gao X, Wan J, Liu B, Ma M, Shen B, Qian SB. Quantitative profiling of initiating ribosomes in vivo. *Nat Methods.* 2015;12(2):147-153. [PMC4344187](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344187/)

32. Lee S, Liu B, Lee S, Huang SX, Shen B, Qian SB. Global mapping of translation initiation sites in mammalian cells at single-nucleotide resolution. *Proc Natl Acad Sci USA.* 2012;109(37):E2424-E2432. [PMC3443142](https://pmc.ncbi.nlm.nih.gov/articles/PMC3443142/)

33. Johnstone TG, Bazzini AA, Giraldez AJ. Upstream ORFs are prevalent translational repressors in vertebrates. *EMBO J.* 2016;35(7):706-723. [PMC4818764](https://pmc.ncbi.nlm.nih.gov/articles/PMC4818764/)

34. Chew GL, Pauli A, Schier AF. Conservation of uORF repressiveness and sequence features in mouse, human, and zebrafish. *Nat Commun.* 2016;7:11663. [PMC4890304](https://pmc.ncbi.nlm.nih.gov/articles/PMC4890304/)

35. Kebaara BW, Atkin AL. Long 3'-UTRs target wild-type mRNAs for nonsense-mediated mRNA decay in Saccharomyces cerevisiae. *Nucleic Acids Res.* 2009;37(9):2771-2778. [PMC2685090](https://pmc.ncbi.nlm.nih.gov/articles/PMC2685090/) *(section re-assigned; see Citation corrections)*

36. Gaba A, Jacobson A, Sachs MS. Ribosome occupancy of the yeast CPA1 upstream open reading frame termination codon modulates nonsense-mediated mRNA decay. *Mol Cell.* 2005;20(5):747-758. [PubMed 16285926](https://pubmed.ncbi.nlm.nih.gov/16285926/)

37. Hinnebusch AG, Ivanov IP, Sonenberg N. Translational control by 5'-untranslated regions of eukaryotic mRNAs. *Science.* 2016;352(6292):1413-1416. [PMC7422601](https://pmc.ncbi.nlm.nih.gov/articles/PMC7422601/)

38. Ferreira JP, Overton KW, Wang CL. Tuning gene expression with synthetic upstream open reading frames. *Proc Natl Acad Sci USA.* 2013;110(28):11284-11289. [PMC3710870](https://pmc.ncbi.nlm.nih.gov/articles/PMC3710870/)

39. Kozak M. Initiation of translation in prokaryotes and eukaryotes. *Gene.* 1999;234(2):187-208. [PubMed 10395892](https://pubmed.ncbi.nlm.nih.gov/10395892/)

40. Lykke-Andersen J, Shu MD, Steitz JA. Communication of the position of exon-exon junctions to the mRNA surveillance machinery by the protein RNPS1. *Science.* 2001;293(5536):1836-1839. [PubMed 11546874](https://pubmed.ncbi.nlm.nih.gov/11546874/)
