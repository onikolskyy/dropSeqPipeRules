import collections
from bin.helperClasses.locus_function import LocusFunctions
from typing import Tuple, List, Dict, Set
from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from bin.helperClasses.gene import Gene
from bin.helperClasses.tags import Tags





def tag_read_with_functional_data(bam_read, gene_interval_tree):
    blocks = bam_read.get_blocks()
    tmp_result = {}
    for b in blocks:
        tmp_result[b] = get_functional_data_for_interval(b, gene_interval_tree)
    final = simplify_functional_data(tmp_result)
    filtered_genes = filter_gene_ids(final.keys(), bam_read, gene_interval_tree)
    filtered_ids = sorted([gene.name for gene in filtered_genes])
    # gene_ids = filter(lambda gene_id: gene_interval_tree.genes[gene_id].is_negative_strand() == bam_read.is_reverse, list(final.keys()))
    # lf_lists = [sorted(list(final[gene_id]), reverse=True) for gene_id in gene_ids]
    #gene_ids_same_strand = filter(lambda gene_id: gene_interval_tree.genes[gene_id].is_negative_strand == read_is_negative_strand, gene_ids)
    #   bam_read.set_tag(Tags.tags_dict["GENE_FUNCTION_TAG"],  ",".join(lf.name for lf_list in lf_lists for lf in lf_list))
    bam_read.set_tag(Tags.tags_dict["GENE_NAME_TAG"], ",".join(filtered_ids))
    #  bam_read.set_tag(Tags.tags_dict["GENE_STRAND_TAG"],   ",".join(["-" if gene_interval_tree.genes[gene_id].is_negative_strand() else "+" for gene_id in gene_ids]))
    #  bam_read.set_tag(Tags.tags_dict["READ_FUNCTION_TAG"], get_best_lf_name(final, gene_ids_same_strand))


def get_functional_data_for_interval(block, gene_interval_tree):
    # type: (Tuple, GeneIntervalTree) -> Dict[str, Set[LocusFunctions]]

    overlap_ids = gene_interval_tree.get_overlaps(block)
    result = {}
    for overlap_id in overlap_ids:
        result[overlap_id] = get_locus_functions_by_interval(gene_interval_tree.genes[overlap_id], block)
    return result


def get_locus_functions_by_interval(g, block):
    # type: (Gene, Tuple) -> Set[LocusFunctions]

    locus_functions = [LocusFunctions.INTERGENIC for i in range(block[1]-block[0]+1)]
    for t in g.transcripts.values():
        t.assign_locus_function_for_range(block[0], locus_functions)
    return set(locus_functions)


def simplify_functional_data(functional_data_map):
    # type: (Dict[Tuple, Dict[str, Set]]) -> Dict[str,set]

    common_gene_ids = set.intersection(*[set(map_for_block.keys()) for block, map_for_block in functional_data_map.items()])
    result = collections.defaultdict(lambda : set())

    for b in functional_data_map.keys():
        for g in common_gene_ids:
            result[g].update(functional_data_map[b][g])
    return result


def get_best_lf_name(lf_map, gene_ids):
    best_lf = LocusFunctions.INTERGENIC
    for gene_id in gene_ids:
        for lf in lf_map[gene_id]:
            if int(best_lf) > int(lf):
                best_lf = lf
                if int(best_lf) >= 5:
                    return best_lf
    return best_lf.name


def getGenesStrictlyOverlapped(blocks, tree):
    blocks_overlap = {}
    for block in blocks:
        overlapped_ids = tree.get_overlaps(block)
        overlapped_genes = [tree.genes[gene_id] for gene_id in overlapped_ids]
        filtered = filter(lambda gene: gene.start <= block[0] and gene.end >= block[1], overlapped_genes)
        blocks_overlap[block] = set(filtered)
    return set.intersection(*[genes for block, genes in blocks_overlap.items()])


def getGenesWitOverlappedUtr(blocks, genes):
    result = set()
    for gene in genes:
        for block in blocks:
            if checkIfUtrOverlapped(block, gene):
                result.add(gene)
    return result


def getGenesWithOverlappedTranscript(blocks, genes):
    result = set()
    for gene in genes:
        for block in blocks:
            if checkIfTranscriptOverlapped(block, gene):
                result.add(gene)
    return result


def getGenesWithOverlappedExon(blocks, genes):
    result = set()
    for gene in genes:
        for block in blocks:
            if checkIfExonOverlapped(block,gene):
                result.add(gene)
    return result


def getGenesWithOverlappedCoding(blocks, genes):
    result = set()
    for gene in genes:
        for block in blocks:
            if checkIfCodingOverlapped(block, gene):
                result.add(gene)
    return result


def checkIfUtrOverlapped(block, gene):
    for t_name, t in gene.transcripts.items():
        utr_overlap1 = checkIfIntervalsOverlap(block, (t.transcription_start, t.coding_start))
        utr_overlap2 = checkIfIntervalsOverlap(block, (t.coding_end, t.transcription_start))
        if utr_overlap1 or utr_overlap2:
            return True
    return False


def checkIfTranscriptOverlapped(block, gene):
    for t_name, t in gene.transcripts.items():
        if checkIfIntervalsOverlap(block, (t.transcription_start, t.transcription_end)):
            return True
    return False

def checkIfExonOverlapped(block, gene):
    for t in gene.transcripts:
        for ex in gene.transcripts[t].exons:
            if  checkIfIntervalsOverlap(ex, block):
                return True
    return False

def checkIfCodingOverlapped(block, gene):
    for transcript_name, t in gene.transcripts.items():
        if checkIfIntervalsOverlap(block, (t.transcription_start, t.transcription_end)):
            return True
    return False


def filterGenesBySameStrand(read, genes):
    same_strand = set()
    for gene in genes:
        if bool(same_strand):
            return set()

    return same_strand


def checkIfIntervalsOverlap(i1,i2):
    if i2[1] >= i1[1]:
        right = i2
        left = i1
    else:
        right = i1
        left = i2
    return right[0] <= left[1]

def filter_gene_ids(gene_ids, read, gi_tree):
    genes = [gi_tree.genes[gene_id] for gene_id in gene_ids]
    # first filter for genes on same strand
    genes_filtered = filter(lambda gene: (gene.strand < 0) == read.is_reverse, genes)
    # retain only those genes where read overlaps an exon
    #return getGenesWithOverlappedExon(read.get_blocks(), genes_filtered)
    return genes_filtered