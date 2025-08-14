       vcf="results/annotated/{sample}.ann.vcf",
        html="results/annotated/{sample}.snpeff.html"
    params:
        db=SNPEFF_DB
    log: "results/logs/snpeff/{sample}.log"
    shell:
        """
        snpEff ann -v {params.db} {input.vcf} > {output.vcf} \
        -stats {output.html} 2> {log}
        """
