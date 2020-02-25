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
    #filtered_genes = filter_gene_ids(final.keys(), bam_read, gene_interval_tree)
    #filtered_ids = sorted([gene.name for gene in filtered_genes])
    filtered_ids = final.keys()
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


def getGenesWithOverlappedExon(read, genes):
    result = set()
    blocks = read.get_blocks()
    for gene in genes:
        for block in blocks:
            if checkIfExonOverlapped(block,gene):
                result.add(gene)
    return result


def checkIfExonOverlapped(block, gene):
    for t in gene.transcripts:
        for ex in gene.transcripts[t].exons:
            if ex[0] <= block[0] and ex[1] >= block[1]:
                return True
    return False


def filter_gene_ids(gene_ids, read, gi_tree):
    genes = [gi_tree.genes[gene_id] for gene_id in gene_ids]
    # first filter for genes on same strand
    genes_filtered = filter(lambda gene: gene.strand == read.is_reverse ,genes)
    # retain only those genes where read overlaps an exon
    return getGenesWithOverlappedExon(read, genes_filtered)