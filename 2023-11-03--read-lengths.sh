mkdir -p readlengths

(for bioproject in $(ls /Users/jeffkaufman/code/mgs-restricted/bioprojects/) ; do
    for rl in $(aws s3 ls s3://nao-restricted/$bioproject/readlengths/ | \
                    awk '{print $NF}' ); do
        if ! [ -e readlengths/$rl ] ; then
            echo s3://nao-restricted/$bioproject/readlengths/$rl
        fi
    done
 done
 for bioproject in $(ls /Users/jeffkaufman/code/mgs-pipeline/bioprojects/) ; do
    for rl in $(aws s3 ls s3://nao-mgs/$bioproject/readlengths/ | \
                    awk '{print $NF}' ); do
        if ! [ -e readlengths/$rl ] ; then
            echo s3://nao-mgs/$bioproject/readlengths/$rl
        fi
    done
done) | xargs -I {} -P 16 aws s3 cp {} readlengths/
