package org.campagnelab.dl.framework.bed;

import it.unimi.dsi.fastutil.objects.*;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.List;

public class BedLoader {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(BedLoader.class);

    /**
     * Read a BED file with three columns: chromosome, start (zero-based), end (zero-based),

     *
     * @param reader reader over a BED File.
     * @return
     * @throws IOException
     */
    public static BEDRecords loadBedFile(final Reader reader) throws IOException {

        BufferedReader br = null;
        BEDRecords records=new BEDRecords();
        try {
            br = new BufferedReader(reader);
            String line;
            //final String header = br.readLine();

            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#")) {
                    final String[] linearray = line.trim().split("\t");
                    if (linearray.length < 3) {
                        LOG.warn("Annotation file, encountered truncated line, ignoring: " + line);
                        continue;
                    }
                    final String chromosome = linearray[0];
                    //           if(!chromosome.equalsIgnoreCase(chroName)) continue;
                    // note start and end are zero-based in the bed format:
                    final int segmentStart = Integer.parseInt(linearray[1]);
                    final int segmentEnd = Integer.parseInt(linearray[2]);
                    final BEDRecord record = new BEDRecord(chromosome, segmentStart, segmentEnd);
                    records.add(record);
                }
            }
        } finally {
            IOUtils.closeQuietly(reader);
        }

        records.sort();
        return records;
    }
}
