from pathway_finder import Finder


def main():
    finder = Finder(accession='NC_000915.1', coreGene="NP_206977.1", hit=5)
    pathways = finder.finder()
    finder.produce(pathways)

    finder = Finder(accession="AAD07482.1")
    pathways = finder.auxFinder()
    finder.produce(pathways)

    finder = Finder(seq="AUGTTTYRRSTVVVVALLISSTUCCYTADQ")
    pathways = finder.seqFinder()
    finder.produce(pathways)


if __name__ == '__main__':
    main()
