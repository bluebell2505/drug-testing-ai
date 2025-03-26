
import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { ArrowLeft, AlertTriangle, CheckCircle, Shield, HelpCircle } from 'lucide-react';
import { Button, Card, CardContent, CardHeader, CardTitle } from '@/components/ui';

interface ToxicityData {
  type: string;
  status: 'Safe' | 'Toxic' | 'Unknown';
  confidence: number;
}

interface DrugInteraction {
  medication: string;
  risk: 'High' | 'Medium' | 'Low';
  mechanism: string;
}

interface TargetPrediction {
  target: string;
  probability: number;
  confidence: 'High' | 'Medium' | 'Low';
}

interface ProjectReportData {
  projectId: number;
  projectName: string;
  executiveSummary: {
    drugLikenessScore: number;
    lipinskiViolations: number;
    toxicityConcerns: number;
    syntheticAccessibility: number;
    bbbPermeability: string;
    overall: string;
  };
  physicochemicalProperties: {
    molecularWeight: string;
    logP: string;
    hBondDonors: string;
    hBondAcceptors: string;
    rotatableBonds: string;
    tpsa: string;
    aromaticRings: string;
    fractionSp3: string;
  };
  drugLikeness: {
    qedScore: string;
    lipinskiViolations: string;
    classification: string;
  };
  admetPredictions: {
    bioavailability: string;
    permeability: string;
    bbbPenetration: string;
    halfLife: string;
    clearance: string;
  };
  toxicityProfile: ToxicityData[];
  drugPropertyProfile: {
    overallScore: number;
    classification: string;
  };
  targetPredictions: TargetPrediction[];
  drugInteractions: DrugInteraction[];
  syntheticAccessibility: {
    score: string;
    factors: string;
  };
}

// Mock data for project reports
const mockReports: ProjectReportData[] = [
  {
    projectId: 1,
    projectName: "Project Alpha",
    executiveSummary: {
      drugLikenessScore: 0.65,
      lipinskiViolations: 0,
      toxicityConcerns: 1,
      syntheticAccessibility: 4.5,
      bbbPermeability: "Medium",
      overall: "Good drug candidate"
    },
    physicochemicalProperties: {
      molecularWeight: "180.2 g/mol",
      logP: "1.8",
      hBondDonors: "1",
      hBondAcceptors: "3",
      rotatableBonds: "5",
      tpsa: "60.0",
      aromaticRings: "1",
      fractionSp3: "0.40"
    },
    drugLikeness: {
      qedScore: "0.65",
      lipinskiViolations: "0",
      classification: "Good"
    },
    admetPredictions: {
      bioavailability: "60%",
      permeability: "High",
      bbbPenetration: "Medium",
      halfLife: "12 hours",
      clearance: "5 mL/min/kg"
    },
    toxicityProfile: [
      { type: "Hepato Toxicity", status: "Safe", confidence: 0.75 },
      { type: "Cardio Toxicity", status: "Safe", confidence: 0.80 },
      { type: "Nephro Toxicity", status: "Safe", confidence: 0.70 },
      { type: "Neuro Toxicity", status: "Toxic", confidence: 0.55 }
    ],
    drugPropertyProfile: {
      overallScore: 0.72,
      classification: "Strong candidate"
    },
    targetPredictions: [
      { target: "EGFR", probability: 0.85, confidence: "High" },
      { target: "HER2", probability: 0.60, confidence: "Medium" },
      { target: "Dopamine Receptors", probability: 0.30, confidence: "Low" }
    ],
    drugInteractions: [
      { medication: "Warfarin", risk: "High", mechanism: "CYP2C9 metabolism competition" },
      { medication: "Aspirin", risk: "Low", mechanism: "Minimal interaction" }
    ],
    syntheticAccessibility: {
      score: "4.5/10 (Moderate)",
      factors: "Complex ring system, multiple rotatable bonds"
    }
  },
  // Additional mock reports for other projects would go here
];

const ProjectReport: React.FC = () => {
  const { projectId } = useParams<{ projectId: string }>();
  const navigate = useNavigate();
  
  const [report, setReport] = useState<ProjectReportData | null>(null);
  
  // Fetch report data
  useEffect(() => {
    if (projectId) {
      const foundReport = mockReports.find(r => r.projectId.toString() === projectId);
      if (foundReport) {
        setReport(foundReport);
      } else {
        // If no report exists for this project ID, create a generic one
        setReport({
          projectId: parseInt(projectId),
          projectName: `Project ${projectId}`,
          executiveSummary: {
            drugLikenessScore: 0.65,
            lipinskiViolations: 0,
            toxicityConcerns: 1,
            syntheticAccessibility: 4.5,
            bbbPermeability: "Medium",
            overall: "Good drug candidate"
          },
          physicochemicalProperties: {
            molecularWeight: "180.2 g/mol",
            logP: "1.8",
            hBondDonors: "1",
            hBondAcceptors: "3",
            rotatableBonds: "5",
            tpsa: "60.0",
            aromaticRings: "1",
            fractionSp3: "0.40"
          },
          drugLikeness: {
            qedScore: "0.65",
            lipinskiViolations: "0",
            classification: "Good"
          },
          admetPredictions: {
            bioavailability: "60%",
            permeability: "High",
            bbbPenetration: "Medium",
            halfLife: "12 hours",
            clearance: "5 mL/min/kg"
          },
          toxicityProfile: [
            { type: "Hepato Toxicity", status: "Safe", confidence: 0.75 },
            { type: "Cardio Toxicity", status: "Safe", confidence: 0.80 },
            { type: "Nephro Toxicity", status: "Safe", confidence: 0.70 },
            { type: "Neuro Toxicity", status: "Toxic", confidence: 0.55 }
          ],
          drugPropertyProfile: {
            overallScore: 0.72,
            classification: "Strong candidate"
          },
          targetPredictions: [
            { target: "EGFR", probability: 0.85, confidence: "High" },
            { target: "HER2", probability: 0.60, confidence: "Medium" },
            { target: "Dopamine Receptors", probability: 0.30, confidence: "Low" }
          ],
          drugInteractions: [
            { medication: "Warfarin", risk: "High", mechanism: "CYP2C9 metabolism competition" },
            { medication: "Aspirin", risk: "Low", mechanism: "Minimal interaction" }
          ],
          syntheticAccessibility: {
            score: "4.5/10 (Moderate)",
            factors: "Complex ring system, multiple rotatable bonds"
          }
        });
      }
    }
  }, [projectId]);
  
  const getToxicityIcon = (status: 'Safe' | 'Toxic' | 'Unknown') => {
    switch (status) {
      case 'Safe':
        return <CheckCircle className="text-green-500" size={18} />;
      case 'Toxic':
        return <AlertTriangle className="text-red-500" size={18} />;
      case 'Unknown':
        return <HelpCircle className="text-gray-500" size={18} />;
    }
  };
  
  const getRiskBadgeClass = (risk: 'High' | 'Medium' | 'Low') => {
    switch (risk) {
      case 'High':
        return "bg-red-500/20 text-red-400";
      case 'Medium':
        return "bg-yellow-500/20 text-yellow-400";
      case 'Low':
        return "bg-green-500/20 text-green-400";
    }
  };
  
  const getConfidenceBadgeClass = (confidence: 'High' | 'Medium' | 'Low') => {
    switch (confidence) {
      case 'High':
        return "bg-green-500/20 text-green-400";
      case 'Medium':
        return "bg-blue-500/20 text-blue-400";
      case 'Low':
        return "bg-gray-500/20 text-gray-400";
    }
  };
  
  if (!report) {
    return (
      <div className="flex flex-col min-h-screen bg-gray-900">
        <Navbar />
        <div className="flex flex-1">
          <Sidebar />
          <main className="flex-1 p-6 flex items-center justify-center">
            <p className="text-white">Loading report...</p>
          </main>
        </div>
        <Footer />
      </div>
    );
  }
  
  return (
    <div className="flex flex-col min-h-screen bg-gray-900">
      <Navbar />
      
      <div className="flex flex-1">
        <Sidebar />
        
        <main className="flex-1 p-6">
          <div className="flex items-center gap-4 mb-6">
            <Button 
              variant="outline" 
              className="text-white" 
              onClick={() => navigate('/projects')}
            >
              <ArrowLeft size={18} className="mr-2" />
              Back to Projects
            </Button>
            <h1 className="text-white text-3xl font-bold">
              AI Drug Report for {report.projectName}
            </h1>
          </div>
          
          {/* Executive Summary */}
          <Card className="bg-gray-800 border-gray-700 mb-6">
            <CardHeader>
              <CardTitle className="text-white text-xl">Executive Summary</CardTitle>
            </CardHeader>
            <CardContent className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
              <div className="bg-gray-700 p-4 rounded-md">
                <h3 className="text-gray-400 text-sm mb-1">Drug-likeness score (QED)</h3>
                <p className="text-white font-medium text-lg">{report.executiveSummary.drugLikenessScore.toFixed(2)}</p>
              </div>
              
              <div className="bg-gray-700 p-4 rounded-md">
                <h3 className="text-gray-400 text-sm mb-1">Lipinski violations</h3>
                <p className="text-white font-medium text-lg">{report.executiveSummary.lipinskiViolations}</p>
              </div>
              
              <div className="bg-gray-700 p-4 rounded-md">
                <h3 className="text-gray-400 text-sm mb-1">Toxicity concerns</h3>
                <p className="text-white font-medium text-lg">{report.executiveSummary.toxicityConcerns}</p>
              </div>
              
              <div className="bg-gray-700 p-4 rounded-md">
                <h3 className="text-gray-400 text-sm mb-1">Synthetic accessibility</h3>
                <p className="text-white font-medium text-lg">{report.executiveSummary.syntheticAccessibility}/10</p>
              </div>
              
              <div className="bg-gray-700 p-4 rounded-md">
                <h3 className="text-gray-400 text-sm mb-1">BBB permeability</h3>
                <p className="text-white font-medium text-lg">{report.executiveSummary.bbbPermeability}</p>
              </div>
              
              <div className="bg-gray-700 p-4 rounded-md">
                <h3 className="text-gray-400 text-sm mb-1">Overall</h3>
                <p className="text-green-400 font-medium text-lg">{report.executiveSummary.overall}</p>
              </div>
            </CardContent>
          </Card>
          
          {/* Physicochemical Properties */}
          <Card className="bg-gray-800 border-gray-700 mb-6">
            <CardHeader>
              <CardTitle className="text-white text-xl">Physicochemical Properties</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="bg-gray-700 rounded-md overflow-hidden">
                <table className="w-full text-left">
                  <thead>
                    <tr className="border-b border-gray-600">
                      <th className="p-4 text-gray-400">Property</th>
                      <th className="p-4 text-gray-400">Value</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr className="border-b border-gray-600">
                      <td className="p-4 text-white">Molecular Weight</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.molecularWeight}</td>
                    </tr>
                    <tr className="border-b border-gray-600">
                      <td className="p-4 text-white">LogP</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.logP}</td>
                    </tr>
                    <tr className="border-b border-gray-600">
                      <td className="p-4 text-white">H-Bond Donors</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.hBondDonors}</td>
                    </tr>
                    <tr className="border-b border-gray-600">
                      <td className="p-4 text-white">H-Bond Acceptors</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.hBondAcceptors}</td>
                    </tr>
                    <tr className="border-b border-gray-600">
                      <td className="p-4 text-white">Rotatable Bonds</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.rotatableBonds}</td>
                    </tr>
                    <tr className="border-b border-gray-600">
                      <td className="p-4 text-white">TPSA</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.tpsa}</td>
                    </tr>
                    <tr className="border-b border-gray-600">
                      <td className="p-4 text-white">Aromatic Rings</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.aromaticRings}</td>
                    </tr>
                    <tr>
                      <td className="p-4 text-white">Fraction SP3</td>
                      <td className="p-4 text-white">{report.physicochemicalProperties.fractionSp3}</td>
                    </tr>
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
          
          {/* Drug Likeness and ADMET */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            {/* Drug Likeness Assessment */}
            <Card className="bg-gray-800 border-gray-700">
              <CardHeader>
                <CardTitle className="text-white text-xl">Drug Likeness Assessment</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="bg-gray-700 p-4 rounded-md mb-3">
                  <h3 className="text-gray-400 text-sm mb-1">QED Score</h3>
                  <p className="text-white font-medium text-lg">{report.drugLikeness.qedScore}</p>
                </div>
                <div className="bg-gray-700 p-4 rounded-md mb-3">
                  <h3 className="text-gray-400 text-sm mb-1">Lipinski Violations</h3>
                  <p className="text-white font-medium text-lg">{report.drugLikeness.lipinskiViolations}</p>
                </div>
                <div className="bg-gray-700 p-4 rounded-md">
                  <h3 className="text-gray-400 text-sm mb-1">Classification</h3>
                  <p className="text-green-400 font-medium text-lg">{report.drugLikeness.classification}</p>
                </div>
              </CardContent>
            </Card>
            
            {/* ADMET Predictions */}
            <Card className="bg-gray-800 border-gray-700">
              <CardHeader>
                <CardTitle className="text-white text-xl">ADMET Predictions</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="bg-gray-700 p-4 rounded-md mb-3">
                  <h3 className="text-gray-400 text-sm mb-1">Bioavailability</h3>
                  <p className="text-white font-medium text-lg">{report.admetPredictions.bioavailability}</p>
                </div>
                <div className="bg-gray-700 p-4 rounded-md mb-3">
                  <h3 className="text-gray-400 text-sm mb-1">Permeability</h3>
                  <p className="text-white font-medium text-lg">{report.admetPredictions.permeability}</p>
                </div>
                <div className="bg-gray-700 p-4 rounded-md mb-3">
                  <h3 className="text-gray-400 text-sm mb-1">Blood-Brain Barrier Penetration</h3>
                  <p className="text-white font-medium text-lg">{report.admetPredictions.bbbPenetration}</p>
                </div>
                <div className="bg-gray-700 p-4 rounded-md mb-3">
                  <h3 className="text-gray-400 text-sm mb-1">Half-Life</h3>
                  <p className="text-white font-medium text-lg">{report.admetPredictions.halfLife}</p>
                </div>
                <div className="bg-gray-700 p-4 rounded-md">
                  <h3 className="text-gray-400 text-sm mb-1">Clearance</h3>
                  <p className="text-white font-medium text-lg">{report.admetPredictions.clearance}</p>
                </div>
              </CardContent>
            </Card>
          </div>
          
          {/* Toxicity Profile */}
          <Card className="bg-gray-800 border-gray-700 mb-6">
            <CardHeader>
              <CardTitle className="text-white text-xl">Toxicity Profile</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="bg-gray-700 rounded-md overflow-hidden">
                <table className="w-full text-left">
                  <thead>
                    <tr className="border-b border-gray-600">
                      <th className="p-4 text-gray-400">Type</th>
                      <th className="p-4 text-gray-400">Status</th>
                      <th className="p-4 text-gray-400">Confidence</th>
                    </tr>
                  </thead>
                  <tbody>
                    {report.toxicityProfile.map((item, index) => (
                      <tr key={index} className={index < report.toxicityProfile.length - 1 ? "border-b border-gray-600" : ""}>
                        <td className="p-4 text-white">{item.type}</td>
                        <td className="p-4 text-white flex items-center">
                          {getToxicityIcon(item.status)}
                          <span className="ml-2">{item.status}</span>
                        </td>
                        <td className="p-4 text-white">{(item.confidence * 100).toFixed(0)}%</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
          
          {/* Drug Property Profile */}
          <Card className="bg-gray-800 border-gray-700 mb-6">
            <CardHeader>
              <CardTitle className="text-white text-xl">Drug Property Profile</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="bg-gray-700 p-5 rounded-md text-center">
                <h3 className="text-gray-400 text-sm mb-2">Overall drug property score (0-1 scale)</h3>
                <div className="relative h-4 bg-gray-600 rounded-full mb-3">
                  <div 
                    className="absolute top-0 left-0 h-full bg-blue-500 rounded-full" 
                    style={{ width: `${report.drugPropertyProfile.overallScore * 100}%` }}
                  ></div>
                </div>
                <p className="text-white font-bold text-2xl">
                  {report.drugPropertyProfile.overallScore.toFixed(2)} 
                  <span className="text-green-400 text-lg ml-2">({report.drugPropertyProfile.classification})</span>
                </p>
              </div>
            </CardContent>
          </Card>
          
          {/* Target Prediction */}
          <Card className="bg-gray-800 border-gray-700 mb-6">
            <CardHeader>
              <CardTitle className="text-white text-xl">Target Prediction</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="bg-gray-700 rounded-md overflow-hidden">
                <table className="w-full text-left">
                  <thead>
                    <tr className="border-b border-gray-600">
                      <th className="p-4 text-gray-400">Target</th>
                      <th className="p-4 text-gray-400">Probability</th>
                      <th className="p-4 text-gray-400">Confidence</th>
                    </tr>
                  </thead>
                  <tbody>
                    {report.targetPredictions.map((item, index) => (
                      <tr key={index} className={index < report.targetPredictions.length - 1 ? "border-b border-gray-600" : ""}>
                        <td className="p-4 text-white">{item.target}</td>
                        <td className="p-4 text-white">{(item.probability * 100).toFixed(0)}%</td>
                        <td className="p-4">
                          <span className={`px-3 py-1 rounded-full text-xs font-medium ${getConfidenceBadgeClass(item.confidence)}`}>
                            {item.confidence}
                          </span>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
          
          {/* Potential Drug-Drug Interactions */}
          <Card className="bg-gray-800 border-gray-700 mb-6">
            <CardHeader>
              <CardTitle className="text-white text-xl">Potential Drug-Drug Interactions</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="bg-gray-700 rounded-md overflow-hidden">
                <table className="w-full text-left">
                  <thead>
                    <tr className="border-b border-gray-600">
                      <th className="p-4 text-gray-400">Medication</th>
                      <th className="p-4 text-gray-400">Risk Level</th>
                      <th className="p-4 text-gray-400">Mechanism</th>
                    </tr>
                  </thead>
                  <tbody>
                    {report.drugInteractions.map((item, index) => (
                      <tr key={index} className={index < report.drugInteractions.length - 1 ? "border-b border-gray-600" : ""}>
                        <td className="p-4 text-white">{item.medication}</td>
                        <td className="p-4">
                          <span className={`px-3 py-1 rounded-full text-xs font-medium ${getRiskBadgeClass(item.risk)}`}>
                            {item.risk}
                          </span>
                        </td>
                        <td className="p-4 text-white">{item.mechanism}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
          
          {/* Synthetic Accessibility */}
          <Card className="bg-gray-800 border-gray-700 mb-6">
            <CardHeader>
              <CardTitle className="text-white text-xl">Synthetic Accessibility</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="bg-gray-700 p-4 rounded-md mb-3">
                <h3 className="text-gray-400 text-sm mb-1">Score</h3>
                <p className="text-white font-medium text-lg">{report.syntheticAccessibility.score}</p>
              </div>
              <div className="bg-gray-700 p-4 rounded-md">
                <h3 className="text-gray-400 text-sm mb-1">Factors affecting synthesis</h3>
                <p className="text-white font-medium">{report.syntheticAccessibility.factors}</p>
              </div>
            </CardContent>
          </Card>
          
          {/* Back to Projects Button */}
          <div className="mb-8">
            <Button 
              variant="outline" 
              className="text-white" 
              onClick={() => navigate('/projects')}
            >
              <ArrowLeft size={18} className="mr-2" />
              Back to Projects
            </Button>
          </div>
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default ProjectReport;
