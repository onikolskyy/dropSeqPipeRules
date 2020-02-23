import collections
from bin.helperClasses.locus_function import LocusFunctions
from typing import Tuple, List, Dict
from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from bin.helperClasses.gene import Gene
from bin.helperClasses.tags import Tags


def tag_read_with_functional_data(bam_read, gene_interval_tree):
    blocks = bam_read.get_blocks()
    tmp_result = {}
    for b in blocks:
        tmp_result[b] = get_functional_data_for_interval(b, gene_interval_tree)
    final = simplify_functional_data(tmp_result)

    gene_ids = list(final.keys())
    gene_ids.sort()
    lf_lists = [sorted(list(final[gene_id]), reverse=True) for gene_id in gene_ids]
    read_is_negative_strand = bam_read.is_reverse
    gene_ids_same_strand = filter(lambda gene_id: gene_interval_tree.genes[gene_id].is_negative_strand == read_is_negative_strand, gene_ids)
    bam_read.set_tag(Tags.tags_dict["GENE_FUNCTION_TAG"],  ", ".join(lf.name for lf_list in lf_lists for lf in lf_list))
    bam_read.set_tag(Tags.tags_dict["GENE_NAME_TAG"], ", ".join(gene_ids))
    bam_read.set_tag(Tags.tags_dict["GENE_STRAND_TAG"],   ", ".join(["-" if gene_interval_tree.genes[gene_id].is_negative_strand() else "+" for gene_id in gene_ids]))
    bam_read.set_tag(Tags.tags_dict["READ_FUNCTION_TAG"], get_best_lf_name(final, gene_ids_same_strand))


def get_functional_data_for_interval(block, gene_interval_tree):
    # type: (Tuple, GeneIntervalTree) -> Dict[str, List[LocusFunctions]]

    overlap_ids = gene_interval_tree.get_overlaps(block)
    result = {}
    for overlap_id in overlap_ids:
        result[overlap_id] = get_locus_functions_by_interval(gene_interval_tree.genes[overlap_id], block)
    return result


def get_locus_functions_by_interval(g, block):
    # type: (Gene, Tuple) -> List[LocusFunctions]

    locus_functions = [LocusFunctions.INTERGENIC for i in range(block[1]-block[0]+1)]
    for t in g.transcripts.values():
        t.assign_locus_function_for_range(block[0], locus_functions)
    return locus_functions


def simplify_functional_data(functional_data_map):
    # type: (Dict[Tuple, Dict[str, list]]) -> Dict[str,set]

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
