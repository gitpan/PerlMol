use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'              => 'Chemistry::InternalCoords',
    'VERSION_FROM'      => 'InternalCoords.pm',
    'PREREQ_PM'         => {
        "Chemistry::Mol" => "0.20",
        "Scalar::Util" => 0,
    },
);
