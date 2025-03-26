import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { useToast } from '@/hooks/use-toast';
import { 
  ArrowLeft, 
  Check, 
  RefreshCw, 
  Save, 
  FileText 
} from 'lucide-react';
import { 
  Button,
  RadioGroup,
  RadioGroupItem,
  Card,
  CardHeader,
  CardContent,
  CardFooter,
  CardTitle,
  CardDescription,
  Checkbox,
  Label,
} from '@/components/ui';

interface Molecule {
  id: number;
  smiles: string;
  properties: {
    solubility: string;
    bindingAffinity: string;
    toxicity: string;
  };
}

interface OptimizedMolecule extends Molecule {
  comparison: {
    solubility: string;
    bindingAffinity: string;
    toxicity: string;
  };
}

interface Project {
  id: number;
  name: string;
  molecules: Molecule[];
}

const MoleculeOptimization: React.FC = () => {
  const { projectId } = useParams<{ projectId: string }>();
  const navigate = useNavigate();
  const { toast } = useToast();
  
  const [project, setProject] = useState<Project | null>(null);
  const [selectedMolecule, setSelectedMolecule] = useState<number | null>(null);
  const [mode, setMode] = useState<'optimize' | 'generate'>('optimize');
  const [optimizationGoals, setOptimizationGoals] = useState<string[]>([]);
  const [desiredInteraction, setDesiredInteraction] = useState<string>('synergistic');
  const [targetProtein, setTargetProtein] = useState<string>('EGFR');
  const [desiredProperties, setDesiredProperties] = useState<string[]>([]);
  const [optimizedMolecules, setOptimizedMolecules] = useState<OptimizedMolecule[]>([]);
  const [isLoading, setIsLoading] = useState<boolean>(false);

  const mockProjects: Project[] = [
    {
      id: 1,
      name: "Project Alpha",
      molecules: [
        { 
          id: 1, 
          smiles: "CC(=O)OC1=CC=CC=C1C(=O)O", 
          properties: { 
            solubility: "Medium", 
            bindingAffinity: "High for EGFR (Score: 8.5)", 
            toxicity: "Low" 
          } 
        },
        { 
          id: 2, 
          smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 
          properties: { 
            solubility: "Low", 
            bindingAffinity: "Low for EGFR (Score: 2.0)", 
            toxicity: "Moderate" 
          } 
        }
      ]
    },
    {
      id: 2,
      name: "Project Beta",
      molecules: [
        { 
          id: 3, 
          smiles: "C1=CC=C2C(=C1)C=CN2", 
          properties: { 
            solubility: "Medium", 
            bindingAffinity: "Medium for HER2 (Score: 6.2)", 
            toxicity: "Low" 
          } 
        }
      ]
    },
    {
      id: 3,
      name: "Project Gamma",
      molecules: []
    },
    {
      id: 4,
      name: "Project Delta",
      molecules: []
    },
    {
      id: 5,
      name: "Project Epsilon",
      molecules: []
    }
  ];

  useEffect(() => {
    if (projectId) {
      const foundProject = mockProjects.find(p => p.id === parseInt(projectId));
      if (foundProject) {
        setProject(foundProject);
      }
    }
  }, [projectId]);

  const handleMoleculeSelect = (moleculeId: number) => {
    setSelectedMolecule(moleculeId === selectedMolecule ? null : moleculeId);
  };

  const handleOptimizationGoalChange = (goal: string) => {
    setOptimizationGoals(prev => 
      prev.includes(goal) 
        ? prev.filter(g => g !== goal) 
        : [...prev, goal]
    );
  };

  const handleDesiredPropertyChange = (property: string) => {
    setDesiredProperties(prev => 
      prev.includes(property) 
        ? prev.filter(p => p !== property) 
        : [...prev, property]
    );
  };

  const handleGenerateOptimizedMolecule = () => {
    if (!selectedMolecule) {
      toast({
        title: "Error",
        description: "Please select a molecule to proceed.",
        variant: "destructive"
      });
      return;
    }

    if (mode === 'optimize' && optimizationGoals.length === 0) {
      toast({
        title: "Error",
        description: "Please select at least one optimization goal.",
        variant: "destructive"
      });
      return;
    }

    if (mode === 'generate' && desiredProperties.length === 0) {
      toast({
        title: "Error",
        description: "Please select at least one desired property.",
        variant: "destructive"
      });
      return;
    }

    setIsLoading(true);

    setTimeout(() => {
      const selectedMoleculeData = project?.molecules.find(m => m.id === selectedMolecule);
      
      if (selectedMoleculeData) {
        const newMolecule: OptimizedMolecule = {
          id: Date.now(),
          smiles: mode === 'optimize' 
            ? "CC(=O)OC1=CC=CC=C1C(=O)OCH3"
            : "CNC(=O)C1=CC=CC=C1NC(=O)C2=CC=CC=C2",
          properties: {
            solubility: mode === 'optimize' && optimizationGoals.includes('solubility') 
              ? "High" 
              : mode === 'generate' && desiredProperties.includes('highSolubility')
                ? "High"
                : selectedMoleculeData.properties.solubility,
            bindingAffinity: mode === 'optimize' && optimizationGoals.includes('bindingAffinity') 
              ? "Very High for EGFR (Score: 9.8)" 
              : mode === 'generate'
                ? `High for ${targetProtein} (Score: 8.9)`
                : selectedMoleculeData.properties.bindingAffinity,
            toxicity: mode === 'optimize' && optimizationGoals.includes('toxicity') 
              ? "Very Low" 
              : mode === 'generate' && desiredProperties.includes('lowToxicity')
                ? "Very Low"
                : selectedMoleculeData.properties.toxicity
          },
          comparison: {
            solubility: mode === 'optimize' && optimizationGoals.includes('solubility')
              ? `Improved from ${selectedMoleculeData.properties.solubility} to High`
              : "No change",
            bindingAffinity: mode === 'optimize' && optimizationGoals.includes('bindingAffinity')
              ? "Improved binding by 15%"
              : "No change",
            toxicity: mode === 'optimize' && optimizationGoals.includes('toxicity')
              ? `Reduced from ${selectedMoleculeData.properties.toxicity} to Very Low`
              : "No change"
          }
        };

        setOptimizedMolecules([newMolecule, ...optimizedMolecules]);
        
        toast({
          title: "Success",
          description: mode === 'optimize' 
            ? "Molecule successfully optimized." 
            : "Complementary molecule successfully generated."
        });
      }
      
      setIsLoading(false);
    }, 2000);
  };

  const handleSaveToProject = (moleculeId: number) => {
    toast({
      title: "Success",
      description: "Molecule saved to project."
    });
  };

  const handleRegenerate = () => {
    handleGenerateOptimizedMolecule();
  };

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
              Molecule Optimization for {project?.name}
            </h1>
          </div>
          
          <div className="mb-8">
            <h2 className="text-white text-xl font-semibold mb-4">Existing Molecules</h2>
            
            {project?.molecules && project.molecules.length > 0 ? (
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                {project.molecules.map(molecule => (
                  <Card 
                    key={molecule.id} 
                    className={`bg-gray-800 border-gray-700 hover:border-blue-400 cursor-pointer transition-all ${
                      selectedMolecule === molecule.id ? 'border-blue-500 ring-2 ring-blue-500' : ''
                    }`}
                    onClick={() => handleMoleculeSelect(molecule.id)}
                  >
                    <CardHeader className="pb-2">
                      <div className="flex justify-between items-start">
                        <CardTitle className="text-white text-lg">Molecule #{molecule.id}</CardTitle>
                        {selectedMolecule === molecule.id && (
                          <div className="bg-blue-500 text-white p-1 rounded-full">
                            <Check size={16} />
                          </div>
                        )}
                      </div>
                      <CardDescription className="text-gray-400 break-all">
                        {molecule.smiles}
                      </CardDescription>
                    </CardHeader>
                    <CardContent>
                      <div className="bg-gray-700 rounded-md p-4 mb-4 h-40 flex items-center justify-center">
                        <p className="text-gray-400 text-center">Molecule visualization placeholder</p>
                      </div>
                      <div className="space-y-2">
                        <p className="text-gray-300"><span className="font-medium">Solubility:</span> {molecule.properties.solubility}</p>
                        <p className="text-gray-300"><span className="font-medium">Binding Affinity:</span> {molecule.properties.bindingAffinity}</p>
                        <p className="text-gray-300"><span className="font-medium">Toxicity:</span> {molecule.properties.toxicity}</p>
                      </div>
                    </CardContent>
                  </Card>
                ))}
              </div>
            ) : (
              <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl text-center">
                <p className="text-gray-400">No molecules found for this project. Generate a new molecule to get started.</p>
              </div>
            )}
          </div>
          
          <div className="mb-8">
            <h2 className="text-white text-xl font-semibold mb-4">Optimization Options</h2>
            
            <Card className="bg-gray-800 border-gray-700">
              <CardContent className="p-6">
                <div className="mb-6">
                  <h3 className="text-white font-medium mb-3">Mode Selection</h3>
                  <RadioGroup 
                    value={mode} 
                    className="flex flex-col space-y-2"
                    onValueChange={(value) => setMode(value as 'optimize' | 'generate')}
                  >
                    <div className="flex items-center space-x-2">
                      <RadioGroupItem value="optimize" id="optimize" className="text-blue-500" />
                      <Label htmlFor="optimize" className="text-white cursor-pointer">Optimize Existing Molecule</Label>
                    </div>
                    <div className="flex items-center space-x-2">
                      <RadioGroupItem value="generate" id="generate" className="text-blue-500" />
                      <Label htmlFor="generate" className="text-white cursor-pointer">Generate Complementary Molecule</Label>
                    </div>
                  </RadioGroup>
                </div>
                
                {mode === 'optimize' ? (
                  <div>
                    <h3 className="text-white font-medium mb-3">Optimization Goals</h3>
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-2 mb-6">
                      <div className="flex items-center space-x-2">
                        <Checkbox 
                          id="solubility" 
                          className="text-blue-500" 
                          checked={optimizationGoals.includes('solubility')}
                          onCheckedChange={() => handleOptimizationGoalChange('solubility')}
                        />
                        <Label htmlFor="solubility" className="text-white cursor-pointer">
                          Improve Solubility (Target: High)
                        </Label>
                      </div>
                      <div className="flex items-center space-x-2">
                        <Checkbox 
                          id="bindingAffinity" 
                          className="text-blue-500"
                          checked={optimizationGoals.includes('bindingAffinity')}
                          onCheckedChange={() => handleOptimizationGoalChange('bindingAffinity')}
                        />
                        <Label htmlFor="bindingAffinity" className="text-white cursor-pointer">
                          Improve Binding Affinity (Target: High)
                        </Label>
                      </div>
                      <div className="flex items-center space-x-2">
                        <Checkbox 
                          id="toxicity" 
                          className="text-blue-500"
                          checked={optimizationGoals.includes('toxicity')}
                          onCheckedChange={() => handleOptimizationGoalChange('toxicity')}
                        />
                        <Label htmlFor="toxicity" className="text-white cursor-pointer">
                          Reduce Toxicity (Target: Low)
                        </Label>
                      </div>
                      <div className="flex items-center space-x-2">
                        <Checkbox 
                          id="bioavailability" 
                          className="text-blue-500"
                          checked={optimizationGoals.includes('bioavailability')}
                          onCheckedChange={() => handleOptimizationGoalChange('bioavailability')}
                        />
                        <Label htmlFor="bioavailability" className="text-white cursor-pointer">
                          Improve Bioavailability (Target: High)
                        </Label>
                      </div>
                    </div>
                    <Button 
                      className="w-full bg-blue-500 hover:bg-blue-600 text-white"
                      onClick={handleGenerateOptimizedMolecule}
                      disabled={isLoading || !selectedMolecule || optimizationGoals.length === 0}
                    >
                      {isLoading ? (
                        <>
                          <RefreshCw size={18} className="mr-2 animate-spin" />
                          Optimizing Molecule...
                        </>
                      ) : (
                        'Optimize Molecule'
                      )}
                    </Button>
                  </div>
                ) : (
                  <div>
                    <h3 className="text-white font-medium mb-3">Desired Interaction</h3>
                    <div className="grid grid-cols-1 gap-2 mb-6">
                      <div className="bg-gray-700 p-3 rounded-md">
                        <RadioGroup 
                          value={desiredInteraction} 
                          className="flex flex-col space-y-2"
                          onValueChange={setDesiredInteraction}
                        >
                          <div className="flex items-center space-x-2">
                            <RadioGroupItem value="synergistic" id="synergistic" className="text-blue-500" />
                            <Label htmlFor="synergistic" className="text-white cursor-pointer">
                              Synergistic Binding
                            </Label>
                          </div>
                          <div className="flex items-center space-x-2">
                            <RadioGroupItem value="allosteric" id="allosteric" className="text-blue-500" />
                            <Label htmlFor="allosteric" className="text-white cursor-pointer">
                              Allosteric Binding
                            </Label>
                          </div>
                          <div className="flex items-center space-x-2">
                            <RadioGroupItem value="competitive" id="competitive" className="text-blue-500" />
                            <Label htmlFor="competitive" className="text-white cursor-pointer">
                              Competitive Binding
                            </Label>
                          </div>
                        </RadioGroup>
                      </div>
                    </div>
                    
                    <h3 className="text-white font-medium mb-3">Target Protein</h3>
                    <div className="bg-gray-700 p-3 rounded-md mb-6">
                      <RadioGroup 
                        value={targetProtein} 
                        className="flex flex-col space-y-2"
                        onValueChange={setTargetProtein}
                      >
                        <div className="flex items-center space-x-2">
                          <RadioGroupItem value="EGFR" id="EGFR" className="text-blue-500" />
                          <Label htmlFor="EGFR" className="text-white cursor-pointer">EGFR</Label>
                        </div>
                        <div className="flex items-center space-x-2">
                          <RadioGroupItem value="HER2" id="HER2" className="text-blue-500" />
                          <Label htmlFor="HER2" className="text-white cursor-pointer">HER2</Label>
                        </div>
                        <div className="flex items-center space-x-2">
                          <RadioGroupItem value="None" id="None" className="text-blue-500" />
                          <Label htmlFor="None" className="text-white cursor-pointer">None</Label>
                        </div>
                      </RadioGroup>
                    </div>
                    
                    <h3 className="text-white font-medium mb-3">Desired Properties</h3>
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-2 mb-6">
                      <div className="flex items-center space-x-2">
                        <Checkbox 
                          id="highSolubility" 
                          className="text-blue-500"
                          checked={desiredProperties.includes('highSolubility')}
                          onCheckedChange={() => handleDesiredPropertyChange('highSolubility')}
                        />
                        <Label htmlFor="highSolubility" className="text-white cursor-pointer">
                          High Solubility
                        </Label>
                      </div>
                      <div className="flex items-center space-x-2">
                        <Checkbox 
                          id="highBindingAffinity" 
                          className="text-blue-500"
                          checked={desiredProperties.includes('highBindingAffinity')}
                          onCheckedChange={() => handleDesiredPropertyChange('highBindingAffinity')}
                        />
                        <Label htmlFor="highBindingAffinity" className="text-white cursor-pointer">
                          High Binding Affinity
                        </Label>
                      </div>
                      <div className="flex items-center space-x-2">
                        <Checkbox 
                          id="lowToxicity" 
                          className="text-blue-500"
                          checked={desiredProperties.includes('lowToxicity')}
                          onCheckedChange={() => handleDesiredPropertyChange('lowToxicity')}
                        />
                        <Label htmlFor="lowToxicity" className="text-white cursor-pointer">
                          Low Toxicity
                        </Label>
                      </div>
                    </div>
                    
                    <Button 
                      className="w-full bg-blue-500 hover:bg-blue-600 text-white"
                      onClick={handleGenerateOptimizedMolecule}
                      disabled={isLoading || !selectedMolecule || desiredProperties.length === 0}
                    >
                      {isLoading ? (
                        <>
                          <RefreshCw size={18} className="mr-2 animate-spin" />
                          Generating Molecule...
                        </>
                      ) : (
                        'Generate Complementary Molecule'
                      )}
                    </Button>
                  </div>
                )}
              </CardContent>
            </Card>
          </div>
          
          <div>
            <h2 className="text-white text-xl font-semibold mb-4">Results</h2>
            
            {optimizedMolecules.length > 0 ? (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {optimizedMolecules.map(molecule => (
                  <Card key={molecule.id} className="bg-gray-800 border-gray-700">
                    <CardHeader>
                      <CardTitle className="text-white">
                        {mode === 'optimize' ? 'Optimized Molecule' : 'Complementary Molecule'}
                      </CardTitle>
                      <CardDescription className="text-gray-400 break-all">
                        {molecule.smiles}
                      </CardDescription>
                    </CardHeader>
                    <CardContent>
                      <div className="bg-gray-700 rounded-md p-4 mb-4 h-40 flex items-center justify-center">
                        <p className="text-gray-400 text-center">Molecule visualization placeholder</p>
                      </div>
                      
                      <div className="space-y-3">
                        <div className="flex flex-col">
                          <p className="text-gray-300 mb-1"><span className="font-medium">Solubility:</span> {molecule.properties.solubility}</p>
                          {molecule.comparison.solubility !== "No change" && (
                            <p className="text-green-400 text-sm">{molecule.comparison.solubility}</p>
                          )}
                        </div>
                        
                        <div className="flex flex-col">
                          <p className="text-gray-300 mb-1"><span className="font-medium">Binding Affinity:</span> {molecule.properties.bindingAffinity}</p>
                          {molecule.comparison.bindingAffinity !== "No change" && (
                            <p className="text-green-400 text-sm">{molecule.comparison.bindingAffinity}</p>
                          )}
                        </div>
                        
                        <div className="flex flex-col">
                          <p className="text-gray-300 mb-1"><span className="font-medium">Toxicity:</span> {molecule.properties.toxicity}</p>
                          {molecule.comparison.toxicity !== "No change" && (
                            <p className="text-green-400 text-sm">{molecule.comparison.toxicity}</p>
                          )}
                        </div>
                      </div>
                    </CardContent>
                    <CardFooter className="flex flex-wrap gap-2">
                      <Button 
                        className="flex-1 gap-1 bg-blue-500 hover:bg-blue-600 text-white"
                        onClick={() => handleSaveToProject(molecule.id)}
                      >
                        <Save size={16} />
                        Save to Project
                      </Button>
                      <Button 
                        variant="outline" 
                        className="flex-1 gap-1 text-white border-gray-600"
                        onClick={() => navigate(`/molecule/${molecule.id}/report`)}
                      >
                        <FileText size={16} />
                        View Report
                      </Button>
                      <Button 
                        variant="outline" 
                        className="flex-1 gap-1 text-white border-gray-600"
                        onClick={handleRegenerate}
                      >
                        <RefreshCw size={16} />
                        Regenerate
                      </Button>
                    </CardFooter>
                  </Card>
                ))}
              </div>
            ) : (
              <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl text-center">
                <p className="text-gray-400">
                  {project?.molecules.length 
                    ? "Select a molecule and optimization options to generate results."
                    : "No molecules available. Add molecules to your project first."}
                </p>
              </div>
            )}
          </div>
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default MoleculeOptimization;
