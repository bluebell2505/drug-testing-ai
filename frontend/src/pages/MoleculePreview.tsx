
import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { ArrowLeft, Settings } from 'lucide-react';
import { Button, Card, CardContent, CardHeader, CardTitle } from '@/components/ui';

const MoleculePreview: React.FC = () => {
  const { projectId } = useParams<{ projectId: string }>();
  const navigate = useNavigate();
  
  // Mock project data
  const [project, setProject] = useState<{ id: number; name: string } | null>(null);
  
  // Mock molecule data for the preview
  const [molecule, setMolecule] = useState<{
    smiles: string;
    molecularWeight: string;
    solubility: string;
    bindingAffinity: string;
    toxicity: string;
    bioavailability: string;
  } | null>(null);
  
  // Mock projects data
  const mockProjects = [
    { id: 1, name: "Project Alpha", molecules: [
      { 
        smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
        molecularWeight: "180.2 g/mol",
        solubility: "Medium",
        bindingAffinity: "High for EGFR (Score: 8.5)",
        toxicity: "Low",
        bioavailability: "High"
      }
    ]},
    { id: 2, name: "Project Beta", molecules: [
      {
        smiles: "C1=CC=C2C(=C1)C=CN2",
        molecularWeight: "150.2 g/mol",
        solubility: "Medium",
        bindingAffinity: "Medium for HER2 (Score: 6.2)",
        toxicity: "Low",
        bioavailability: "High"
      }
    ]},
    { id: 3, name: "Project Gamma", molecules: [] },
    { id: 4, name: "Project Delta", molecules: [] },
    { id: 5, name: "Project Epsilon", molecules: [] }
  ];
  
  // Fetch project data
  useEffect(() => {
    if (projectId) {
      const foundProject = mockProjects.find(p => p.id.toString() === projectId);
      if (foundProject) {
        setProject(foundProject);
        // Set the first molecule of the project for preview (if any)
        if (foundProject.molecules && foundProject.molecules.length > 0) {
          setMolecule(foundProject.molecules[0]);
        }
      }
    }
  }, [projectId]);
  
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
              Molecule Preview for {project?.name}
            </h1>
          </div>
          
          {molecule ? (
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              {/* 3D Visualization Section */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader>
                  <CardTitle className="text-white">3D Visualization</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="bg-gray-700 rounded-md p-6 h-80 flex items-center justify-center">
                    <div className="text-center">
                      <p className="text-gray-400 mb-3">3D Visualization of {molecule.smiles}</p>
                      <p className="text-gray-500 text-sm">
                        Placeholder for 3D molecular visualization <br />
                        (A 3D model would be rendered here in a production environment)
                      </p>
                    </div>
                  </div>
                </CardContent>
              </Card>
              
              {/* Key Characteristics Section */}
              <Card className="bg-gray-800 border-gray-700">
                <CardHeader>
                  <CardTitle className="text-white">Key Characteristics</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="space-y-4">
                    <div className="bg-gray-700 p-4 rounded-md">
                      <h3 className="text-white font-medium mb-1">Molecular Weight</h3>
                      <p className="text-blue-400 text-lg font-semibold">{molecule.molecularWeight}</p>
                    </div>
                    
                    <div className="bg-gray-700 p-4 rounded-md">
                      <h3 className="text-white font-medium mb-1">Solubility</h3>
                      <p className="text-blue-400 text-lg font-semibold">{molecule.solubility}</p>
                    </div>
                    
                    <div className="bg-gray-700 p-4 rounded-md">
                      <h3 className="text-white font-medium mb-1">Binding Affinity</h3>
                      <p className="text-blue-400 text-lg font-semibold">{molecule.bindingAffinity}</p>
                    </div>
                    
                    <div className="bg-gray-700 p-4 rounded-md">
                      <h3 className="text-white font-medium mb-1">Toxicity</h3>
                      <p className="text-blue-400 text-lg font-semibold">{molecule.toxicity}</p>
                    </div>
                    
                    <div className="bg-gray-700 p-4 rounded-md">
                      <h3 className="text-white font-medium mb-1">Bioavailability</h3>
                      <p className="text-blue-400 text-lg font-semibold">{molecule.bioavailability}</p>
                    </div>
                  </div>
                </CardContent>
              </Card>
              
              {/* Action Button */}
              <div className="lg:col-span-2">
                <Button 
                  className="w-full mt-4 bg-blue-500 hover:bg-blue-600 text-white p-6 text-lg"
                  onClick={() => navigate(`/projects/${projectId}/molecule-optimization`)}
                >
                  <Settings size={20} className="mr-2" />
                  Optimize This Molecule
                </Button>
              </div>
            </div>
          ) : (
            <div className="bg-gray-800/50 border border-gray-700 p-8 rounded-xl text-center">
              <p className="text-gray-400 mb-6">No molecules found for this project. Generate a new molecule to get started.</p>
              <Button 
                className="bg-blue-500 hover:bg-blue-600 text-white"
                onClick={() => navigate(`/projects/${projectId}/molecule-optimization`)}
              >
                Generate New Molecule
              </Button>
            </div>
          )}
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default MoleculePreview;
