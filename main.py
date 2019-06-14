from pathway_finder import Finder


def main():
    finder = Finder(accession='AL123456.3', coreGene="Rv0032")
    pathways = finder.finder()
    finder.produce(pathways)


if __name__ == '__main__':
    main()
