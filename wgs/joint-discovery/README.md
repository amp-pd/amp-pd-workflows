The workflow in this repository is the Joint Discovery workflow provided
by the Broad Institute, along with a handful of small modifications to allow
the workflow to run more reliably for cohorts of around 4,000 WGS samples.

A copy of the joint-discovery-gatk4.wdl file was made on July 1, 2019
from the Broad's github repo:

https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/240752418d7a1f64aaaf73ff425e22a45969519f/joint-discovery-gatk4.wdl

and imported into the FireCloud Methods Repository.


A few updates needed to be made to get it running reliably in Terra for AMP PD:

1- Allow for passing a "user_project_id" such that ImportGVCFs can support
requester pays buckets.

2- Add support for maxRetries, as shards can arbitrarily fail due to
Pipelines API Error 10.
See https://support.terra.bio/hc/en-us/community/posts/360046714292.

3- Force the number of shards to be under 10,000 such that the workflow
doesn't create more than 50,000 nodes.
See https://support.terra.bio/hc/en-us/community/posts/360048098952.

Setting the `merge_count` to 3 results in 6,798 shards.
Without explicitly setting it, the merge_count is 2, resulting in 10,187.

4- Explicitly use a full-priced VM for SNPsVariantRecalibrator.
This is a single gating node for the next batch of sharded steps.
Prefer not to use preemptible VMs here.

5- Increase memory available for VariantRecalibrator.
For runs of around 4,000 samples, the 3 GB default Java memory was nowhere
near sufficient. With higher values, like 32 GB, this step finished reliably
and much more quickly (under 4 hrs).

6- Correct the path to gsutil for metrics gathering.
See https://github.com/gatk-workflows/gatk4-germline-snps-indels/issues/41.
