package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

import java.util.function.Function;

public interface FillInFeaturesFunction extends Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> {

}
