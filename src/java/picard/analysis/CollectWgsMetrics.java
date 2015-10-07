package picard.analysis;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SamLocusIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.MathUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.
 *
 * @author tfennell
 */
@CommandLineProgramProperties(
        usage = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
                "whole genome sequencing experiments.",
        usageShort = "Writes whole genome sequencing-related metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectWgsMetrics extends CommandLineProgram {

    public static final int ARRAYLENGTH = 1000;
    private static final byte BAD_QUALITY = -1;
    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output metrics file.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "The reference sequence fasta aligned to.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName = "MQ", doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.", overridable = true)
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(shortName = "CAP", doc = "Treat bases with coverage exceeding this value as if they had coverage at this value.", overridable = true)
    public int COVERAGE_CAP = 250;

    @Option(doc = "For debugging purposes, stop after processing this many genomic bases.")
    public long STOP_AFTER = -1;

    @Option(doc = "Determines whether to include the base quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

    @Option(doc = "If true, count unpaired reads, and paired reads with one end unmapped")
    public boolean COUNT_UNPAIRED = false;

    private final Log log = Log.getInstance(CollectWgsMetrics.class);

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class WgsMetrics extends MetricBase {
        /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
        public long GENOME_TERRITORY;
        /** The mean coverage in bases of the genome territory, after all filters are applied. */
        public double MEAN_COVERAGE;
        /** The standard deviation of coverage of the genome after all filters are applied. */
        public double SD_COVERAGE;
        /** The median coverage in bases of the genome territory, after all filters are applied. */
        public double MEDIAN_COVERAGE;
        /** The median absolute deviation of coverage of the genome after all filters are applied. */
        public double MAD_COVERAGE;

        /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
        public double PCT_EXC_MAPQ;
        /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
        public double PCT_EXC_DUPE;
        /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
        public double PCT_EXC_UNPAIRED;
        /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
        public double PCT_EXC_BASEQ;
        /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
        public double PCT_EXC_OVERLAP;
        /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
        public double PCT_EXC_CAPPED;
        /** The total fraction of aligned bases excluded due to all filters. */
        public double PCT_EXC_TOTAL;

        /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
        public double PCT_5X;
        /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
        public double PCT_10X;
        /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
        public double PCT_15X;
        /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
        public double PCT_20X;
        /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
        public double PCT_25X;
        /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
        public double PCT_30X;
        /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
        public double PCT_40X;
        /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
        public double PCT_50X;
        /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
        public double PCT_60X;
        /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
        public double PCT_70X;
        /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
        public double PCT_80X;
        /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
        public double PCT_90X;
        /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
        public double PCT_100X;
    }

    public static void main(final String[] args) {
        new CollectWgsMetrics().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        // Setup all the inputs
        final ProgressLogger progress = new ProgressLogger(log, 10000000, "Processed", "loci");
        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SamLocusIterator iterator = getLocusIterator(in);

        final List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
        final CountingFilter dupeFilter = new CountingDuplicateFilter();
        final CountingFilter mapqFilter = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
        final CountingPairedFilter pairFilter = new CountingPairedFilter();
        filters.add(mapqFilter);
        filters.add(dupeFilter);
        if (!COUNT_UNPAIRED) {
            filters.add(pairFilter);
        }
        filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
        iterator.setSamFilters(filters);
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
        iterator.setQualityScoreCutoff(0);        // Handled separately because we want to count bases
        iterator.setIncludeNonPfReads(false);

        final int max = COVERAGE_CAP;
        final long[] HistogramArray = new long[max + 1];
        long[] baseQHistogramArray = new long[Byte.MAX_VALUE];
        final boolean usingStopAfter = STOP_AFTER > 0;
        final long stopAfter = STOP_AFTER - 1;
        long counter = 0;

        long basesExcludedByBaseq = 0;
        long basesExcludedByOverlap = 0;
        long basesExcludedByCapping = 0;

        class CWgsQualities {
            private int[] _countBasesExcludedByBaseq;
            private int[] _countBasesExcludedByOverlap;

            private int[] _passCount;
            private byte[][] _qualities;
            private long[] _qualitiesCount;
            private int _length;
            int _index;
            int _pass;
            private long[] _tmpQualitiesCount;
            private int[] _readNamesSize;

            public CWgsQualities(int length) {
                _countBasesExcludedByBaseq = new int[length];
                _countBasesExcludedByOverlap = new int[length];
                _passCount = new int[length];
                _qualities = new byte[length][];
                _readNamesSize = new  int[length];
                _tmpQualitiesCount = new long[Byte.MAX_VALUE];
                _qualitiesCount = new long[Byte.MAX_VALUE];
                _length = length;

            }

            public int getIndex(int i) {
                _index = i % _length;
                _pass = i / _length;
                if (_passCount[_index] != _pass) {
                    _countBasesExcludedByBaseq[_index] = 0;
                    _countBasesExcludedByOverlap[_index] = 0;
                    _readNamesSize[_index] = 0;
                    _passCount[_index]=_pass;

                }

                return _index;
            }

            public int getCountBasesExcludedByBaseq(int i) {
//                int index = i % _length;
//                int pass = i / _length;
                if (_passCount[_index] != _pass) {
                    System.out.println("Error");
                }
                return _countBasesExcludedByBaseq[i];
            }

            public void incrimentCountBasesExcludedByBaseq(int i) {

//                int pass = i / _length;
//                int index = i % _length;

//                if (_passCount[index] != pass) {
//                    _passCount[index] = pass;
//                    _countBasesExcludedByBaseq[index] = 0;
//                }
                _countBasesExcludedByBaseq[i]++;

            }

            public int getCountBasesExcludedByOverlap(int i) {
//                int index = i % _length;
//                int pass = i / _length;
                if (_passCount[_index] != _pass) {
                    System.out.println("Error");
                }
                return _countBasesExcludedByOverlap[i];
            }

            public void incrimentCountBasesExcludedByOverlap(int i) {

//                int pass = i / _length;
//                int index = i % _length;
//
//                if (_passCount[index] != pass) {
//                    _passCount[index] = pass;
//                    _countBasesExcludedByBaseq[index] = 0;
//                }
                _countBasesExcludedByOverlap[i]++;

            }


            public byte[][] getQualities() {
                return _qualities;
            }

            public void setQualitiesForOne(int offset, int i, byte value, final int length) {
                if (_qualities[offset] == null) _qualities[offset] = new byte[length];
                _qualities[offset][i] = value;
            }

            public void setTmpQualities(int i) {
                _tmpQualitiesCount[i]++;


            }

            public void clearTmpQualities() {
                _tmpQualitiesCount = new long[Byte.MAX_VALUE];
            }

            public void calculateQualities() {
                for (int i = 0; i < _tmpQualitiesCount.length; i++) {
                    _qualitiesCount[i]+= _tmpQualitiesCount[i];
                }
                clearTmpQualities();
            }


//            public void setQualities(int offset, byte[] qualities, int length) {
//                _qualities[offset] = new byte[length];
//                for (int i = 0; i < _qualities.length; i++) {
//                    _qualities[offset][i] = qualities[i];
//                }
//            }

            public int length() {
                return _length;
            }

            public long[] getQualitiesCount() {
                return _qualitiesCount;
            }

            public void incrimentReadNameSize(final int i) {
                _readNamesSize[i]++;
            }

            public int getReadNamesSize(final int i) {
                return  _readNamesSize[i];
            }


//            public void shiftArray(final int position) {
//
//                List<Integer> lst = new ArrayList<>(_length);
//
//
//                int tmpCount[] = new int[_length];
//                byte[][] tmpQuality = new byte[_length][];
////                System.arraycopy(_countBasesExcludedByBaseq, position, tmpCount, 0, _length - position);
////                System.arraycopy(_qualities, position, tmpQuality, 0, _length - position);
//                for (int i = position; i < _length; i++) {
//                    tmpCount[i - position] = _countBasesExcludedByBaseq[i];
//                    tmpQuality[i - position] = _qualities[i];
//                }
//                _countBasesExcludedByBaseq = tmpCount;
//                _qualities = tmpQuality;
//            }
        }

//        int shift = 0;
        CWgsQualities cWgsQualities = new CWgsQualities(ARRAYLENGTH);
        // Loop through all the loci
        while (iterator.hasNext()) {

            final SamLocusIterator.LocusInfo info = iterator.next();

            // Check that the reference is not N
            final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (base == 'N') continue;


            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<String>(info.getRecordAndPositions().size());
            int pileupSize = 0;

            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {

                if (!recs.isRecordProcessed()) {

                    final int length = recs.getRecord().getBaseQualities().length;

                    boolean pileupsSizeIncrmentIsFirst = true;
                    boolean basesExcludedByOverlapIsFirst = true;
                    boolean basesIsExcludedByOverlap = false;


                    for (int i = recs.getOffset(); i < length; i++) {
                        final byte quality = recs.getRecord().getBaseQualities()[i];
                        int index = cWgsQualities.getIndex(i-recs.getOffset()+info.getPosition());
                        if (quality < MINIMUM_BASE_QUALITY) {
//                            cWgsQualities.incrimentCountBasesExcludedByBaseq(i - recs.getOffset() + info.getPosition() - shift);
                            cWgsQualities.incrimentCountBasesExcludedByBaseq(index);
//  cWgsQualities.setQualitiesForOne(info.getPosition() - shift, i, BAD_QUALITY, length);
                        } else {

                            if (basesExcludedByOverlapIsFirst) {
                                basesIsExcludedByOverlap = readNames.add(recs.getRecord().getReadName());
                                basesExcludedByOverlapIsFirst = false;
                            }
                            if (!basesIsExcludedByOverlap) {
                                cWgsQualities.incrimentCountBasesExcludedByOverlap(index);
                            } else {
                                cWgsQualities.incrimentReadNameSize(index);
                                if (pileupsSizeIncrmentIsFirst) {
                                    pileupSize++;
                                    pileupsSizeIncrmentIsFirst = false;
                                }
                                if (pileupSize<= max) {
                                    cWgsQualities.setTmpQualities(quality);
                                }

                            }

                            //TODO 1 raz +
//                            if (pileupsSizeIncrmentIsFirst) {
//                            }
//                            if (pileupSize <= max) {
//
//                                final byte[] qualities = cWgsQualities.getQualities()[info.getPosition() - shift];
//                                if (qualities != null) {
//                                    if (qualities[recs.getOffset()] != BAD_QUALITY) baseQHistogramArray[qualities[recs.getOffset()]]++;
//                                }
//                            }
//TODO add Array overlap
//                            tmpbaseQHistogramArray[quality]++;

//                            cWgsQualities.setQualitiesForOne(info.getPosition() - shift, i, quality, length);
                        }
                    }
                    recs.processRecord();
                }

                /*if (cWgsQualities.getCountBasesExcludedByBaseq()[info.getPosition() + recs.getOffset() - shift] > 0) {
                    basesExcludedByCapping += cWgsQualities.getCountBasesExcludedByBaseq()[info.getPosition() + recs.getOffset() - shift];
                    continue;
                }*/

//                if (!readNames.add(recs.getRecord().getReadName())) {
//                    ++basesExcludedByOverlap;//TODO iz massiva
////                    continue;
//                }
            }

//            if (cWgsQualities.getCountBasesExcludedByBaseq(info.getPosition()) > 0) {
//                basesExcludedByCapping += cWgsQualities.getCountBasesExcludedByBaseq(info.getPosition());
////                continue;
//            }
            int index = cWgsQualities.getIndex(info.getPosition());
            basesExcludedByBaseq += cWgsQualities.getCountBasesExcludedByBaseq(index);
            basesExcludedByOverlap += cWgsQualities.getCountBasesExcludedByOverlap(index);
//            cWgsQualities.calculateQualities();


            final int readNamesSize = cWgsQualities.getReadNamesSize(index);
            final int depth = Math.min(readNamesSize, max);
            if (depth < readNamesSize) basesExcludedByCapping += readNamesSize - max;
            HistogramArray[depth]++;

            // Record progress and perhaps stop
            progress.record(info.getSequenceName(), info.getPosition());
            if (usingStopAfter && ++counter > stopAfter) break;
        }
        baseQHistogramArray = cWgsQualities.getQualitiesCount();

//        System.out.println("basesExcludedByBaseq: "+basesExcludedByBaseq);
//        System.out.println("basesExcludedByOverlap: "+basesExcludedByOverlap);
////        System.out.println("depth: "+ depth);
//        System.out.println("baseQHistogramArray:");
//        for (int i = 0;i<baseQHistogramArray.length;i++){
//            System.out.print(baseQHistogramArray[i] + " ");
//        }
//        System.out.println();
//
//        System.out.println("HistogramArray:");
//
//        for (int i = 0;i<HistogramArray.length;i++){
//            System.out.print(HistogramArray[i]+" ");
//        }


        // Construct and write the outputs
        final Histogram<Integer> histo = new Histogram<Integer>("coverage", "count");
        for (int i = 0; i < HistogramArray.length; ++i) {
            histo.increment(i, HistogramArray[i]);
        }

        // Construct and write the outputs
        final Histogram<Integer> baseQHisto = new Histogram<Integer>("value", "baseq_count");
        for (int i = 0; i < baseQHistogramArray.length; ++i) {
            baseQHisto.increment(i, baseQHistogramArray[i]);
        }

        final WgsMetrics metrics = generateWgsMetrics();
        metrics.GENOME_TERRITORY = (long) histo.getSumOfValues();
        metrics.MEAN_COVERAGE = histo.getMean();
        metrics.SD_COVERAGE = histo.getStandardDeviation();
        metrics.MEDIAN_COVERAGE = histo.getMedian();
        metrics.MAD_COVERAGE = histo.getMedianAbsoluteDeviation();

        final long basesExcludedByDupes = getBasesExcludedBy(dupeFilter);
        final long basesExcludedByMapq = getBasesExcludedBy(mapqFilter);
        final long basesExcludedByPairing = getBasesExcludedBy(pairFilter);
        final double total = histo.getSum();
        final double totalWithExcludes = total + basesExcludedByDupes + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq + basesExcludedByOverlap + basesExcludedByCapping;
        metrics.PCT_EXC_DUPE = basesExcludedByDupes / totalWithExcludes;
        metrics.PCT_EXC_MAPQ = basesExcludedByMapq / totalWithExcludes;
        metrics.PCT_EXC_UNPAIRED = basesExcludedByPairing / totalWithExcludes;
        metrics.PCT_EXC_BASEQ = basesExcludedByBaseq / totalWithExcludes;
        metrics.PCT_EXC_OVERLAP = basesExcludedByOverlap / totalWithExcludes;
        metrics.PCT_EXC_CAPPED = basesExcludedByCapping / totalWithExcludes;
        metrics.PCT_EXC_TOTAL = (totalWithExcludes - total) / totalWithExcludes;

        metrics.PCT_5X = MathUtil.sum(HistogramArray, 5, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_10X = MathUtil.sum(HistogramArray, 10, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_15X = MathUtil.sum(HistogramArray, 15, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_20X = MathUtil.sum(HistogramArray, 20, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_25X = MathUtil.sum(HistogramArray, 25, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_30X = MathUtil.sum(HistogramArray, 30, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_40X = MathUtil.sum(HistogramArray, 40, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_50X = MathUtil.sum(HistogramArray, 50, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_60X = MathUtil.sum(HistogramArray, 60, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_70X = MathUtil.sum(HistogramArray, 70, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_80X = MathUtil.sum(HistogramArray, 80, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_90X = MathUtil.sum(HistogramArray, 90, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;
        metrics.PCT_100X = MathUtil.sum(HistogramArray, 100, HistogramArray.length) / (double) metrics.GENOME_TERRITORY;

        final MetricsFile<WgsMetrics, Integer> out = getMetricsFile();
        out.addMetric(metrics);
        out.addHistogram(histo);
        if (INCLUDE_BQ_HISTOGRAM) {
            out.addHistogram(baseQHisto);
        }
        out.write(OUTPUT);

        return 0;
    }

    protected WgsMetrics generateWgsMetrics() {
        return new WgsMetrics();
    }

    protected long getBasesExcludedBy(final CountingFilter filter) {
        return filter.getFilteredBases();
    }

    protected SamLocusIterator getLocusIterator(final SamReader in) {
        return new SamLocusIterator(in);
    }
}

/**
 * A SamRecordFilter that counts the number of aligned bases in the reads which it filters out. Abstract and designed
 * to be subclassed to implement the desired filter.
 */
abstract class CountingFilter implements SamRecordFilter {
    private long filteredRecords = 0;
    private long filteredBases = 0;

    /** Gets the number of records that have been filtered out thus far. */
    public long getFilteredRecords() { return this.filteredRecords; }

    /** Gets the number of bases that have been filtered out thus far. */
    public long getFilteredBases() { return this.filteredBases; }

    @Override
    public final boolean filterOut(final SAMRecord record) {
        final boolean filteredOut = reallyFilterOut(record);
        if (filteredOut) {
            ++filteredRecords;
            for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                this.filteredBases += block.getLength();
            }
        }
        return filteredOut;
    }

    abstract public boolean reallyFilterOut(final SAMRecord record);

    @Override
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        throw new UnsupportedOperationException();
    }
}

/** Counting filter that discards reads that have been marked as duplicates. */
class CountingDuplicateFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getDuplicateReadFlag(); }
}

/** Counting filter that discards reads below a configurable mapping quality threshold. */
class CountingMapQFilter extends CountingFilter {
    private final int minMapq;

    CountingMapQFilter(final int minMapq) { this.minMapq = minMapq; }

    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getMappingQuality() < minMapq; }
}

/** Counting filter that discards reads that are unpaired in sequencing and paired reads who's mates are not mapped. */
class CountingPairedFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return !record.getReadPairedFlag() || record.getMateUnmappedFlag(); }
}

