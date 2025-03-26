
import React, { useState } from 'react';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { Button } from '@/components/ui/button';
import { Beaker, Atom, GitBranch, Database, ArrowRight } from 'lucide-react';
import { useToast } from "@/hooks/use-toast";
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { useIsMobile } from '@/hooks/use-mobile';

const MoleculeGenerator: React.FC = () => {
  const [isGenerating, setIsGenerating] = useState(false);
  const [smiles, setSmiles] = useState('');
  const [targetProtein, setTargetProtein] = useState('');
  const [solubility, setSolubility] = useState('');
  const { toast } = useToast();
  const isMobile = useIsMobile();

  // Advanced options state
  const [admetProperties, setAdmetProperties] = useState(false);
  const [bindingAffinities, setBindingAffinities] = useState(false);
  const [deepLearning, setDeepLearning] = useState(false);
  const [molecularFingerprints, setMolecularFingerprints] = useState(false);
  const [physicalProperties, setPhysicalProperties] = useState(false);

  const handleGenerate = (e: React.FormEvent) => {
    e.preventDefault();
    setIsGenerating(true);
    
    // Simulate API call
    setTimeout(() => {
      toast({
        title: "Success",
        description: "Molecule generated successfully!",
      });
      setIsGenerating(false);
    }, 2000);
  };

  return (
    <div className="flex flex-col min-h-screen bg-gray-900">
      <Navbar />
      
      <div className="flex flex-1">
        <Sidebar />
        
        <main className="flex-1 p-6 overflow-auto">
          <h1 className="text-white text-3xl font-bold mb-6">Molecule Generator</h1>
          
          <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl mb-6">
            <form onSubmit={handleGenerate}>
              <div className="mb-4">
                <Label htmlFor="smiles" className="block text-white font-medium mb-2">Enter SMILES or Name</Label>
                <Input 
                  id="smiles" 
                  value={smiles}
                  onChange={(e) => setSmiles(e.target.value)}
                  className="w-full bg-gray-700 border border-gray-600 text-white p-3 rounded-lg"
                  placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O"
                />
                <p className="text-gray-400 text-sm mt-1">
                  AI will extract and analyze information about atom types, bond types, and stereochemistry
                </p>
              </div>
              
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
                <div>
                  <Label htmlFor="target" className="block text-white font-medium mb-2">Target Protein</Label>
                  <select 
                    id="target" 
                    value={targetProtein}
                    onChange={(e) => setTargetProtein(e.target.value)}
                    className="w-full bg-gray-700 border border-gray-600 text-white p-3 rounded-lg"
                  >
                    <option value="">None</option>
                    <option value="EGFR">EGFR</option>
                    <option value="HER2">HER2</option>
                    <option value="CDK">CDK</option>
                    <option value="ACE2">ACE2</option>
                  </select>
                </div>
                <div>
                  <Label htmlFor="solubility" className="block text-white font-medium mb-2">Desired Solubility</Label>
                  <select 
                    id="solubility" 
                    value={solubility}
                    onChange={(e) => setSolubility(e.target.value)}
                    className="w-full bg-gray-700 border border-gray-600 text-white p-3 rounded-lg"
                  >
                    <option value="">None</option>
                    <option value="high">High</option>
                    <option value="medium">Medium</option>
                    <option value="low">Low</option>
                  </select>
                </div>
              </div>

              <div className="bg-gray-700/30 border border-gray-600 p-4 rounded-lg mb-6">
                <h3 className="text-white font-medium mb-3">Advanced Properties</h3>
                
                <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 gap-3">
                  <div className="flex items-center">
                    <input 
                      type="checkbox" 
                      id="physical-properties" 
                      className="w-4 h-4 mr-2"
                      checked={physicalProperties}
                      onChange={(e) => setPhysicalProperties(e.target.checked)}
                    />
                    <Label htmlFor="physical-properties" className="text-white text-sm">
                      Physical Properties
                    </Label>
                  </div>
                  <div className="flex items-center">
                    <input 
                      type="checkbox" 
                      id="molecular-fingerprints" 
                      className="w-4 h-4 mr-2"
                      checked={molecularFingerprints}
                      onChange={(e) => setMolecularFingerprints(e.target.checked)}
                    />
                    <Label htmlFor="molecular-fingerprints" className="text-white text-sm">
                      Molecular Fingerprints
                    </Label>
                  </div>
                  <div className="flex items-center">
                    <input 
                      type="checkbox" 
                      id="admet-properties" 
                      className="w-4 h-4 mr-2"
                      checked={admetProperties}
                      onChange={(e) => setAdmetProperties(e.target.checked)}
                    />
                    <Label htmlFor="admet-properties" className="text-white text-sm">
                      ADMET Properties
                    </Label>
                  </div>
                  <div className="flex items-center">
                    <input 
                      type="checkbox" 
                      id="deep-learning" 
                      className="w-4 h-4 mr-2"
                      checked={deepLearning}
                      onChange={(e) => setDeepLearning(e.target.checked)}
                    />
                    <Label htmlFor="deep-learning" className="text-white text-sm">
                      Deep Learning Approaches
                    </Label>
                  </div>
                  <div className="flex items-center">
                    <input 
                      type="checkbox" 
                      id="binding-affinities" 
                      className="w-4 h-4 mr-2"
                      checked={bindingAffinities}
                      onChange={(e) => setBindingAffinities(e.target.checked)}
                    />
                    <Label htmlFor="binding-affinities" className="text-white text-sm">
                      Binding Affinities (TC-DTA)
                    </Label>
                  </div>
                </div>

                {physicalProperties && (
                  <div className="mt-3 pl-4 border-l-2 border-blue-500">
                    <p className="text-gray-300 text-sm">
                      Physical and chemical test attributes:
                    </p>
                    <ul className="text-gray-400 text-xs mt-1 list-disc pl-4 space-y-1">
                      <li>Solubility in different media and solvents</li>
                      <li>Dissolution characteristics of active pharmaceutical ingredients</li>
                      <li>Solid state properties like polymorphism, particle size and shape</li>
                      <li>Partition coefficient</li>
                      <li>pKa determination</li>
                    </ul>
                  </div>
                )}

                {molecularFingerprints && (
                  <div className="mt-3 pl-4 border-l-2 border-blue-500">
                    <p className="text-gray-300 text-sm">
                      Generate appropriate molecular fingerprints such as Morgan fingerprint and MACCS fingerprint,
                      which are essential for encoding chemical structures and predicting molecular properties.
                    </p>
                  </div>
                )}

                {admetProperties && (
                  <div className="mt-3 pl-4 border-l-2 border-blue-500">
                    <p className="text-gray-300 text-sm">
                      ADMET Properties: Absorption, Distribution, Metabolism, Excretion, Toxicity
                    </p>
                  </div>
                )}

                {deepLearning && (
                  <div className="mt-3 pl-4 border-l-2 border-blue-500">
                    <p className="text-gray-300 text-sm">
                      Deep learning approaches, particularly using transformer architectures, 
                      have shown promise in learning molecular representations for property prediction.
                    </p>
                  </div>
                )}

                {bindingAffinities && (
                  <div className="mt-3 pl-4 border-l-2 border-blue-500">
                    <p className="text-gray-300 text-sm">
                      Predict binding affinities between the drug candidate and relevant biological targets
                      using methodology like the TC-DTA model, which leverages convolutional neural networks
                      and transformer architecture.
                    </p>
                  </div>
                )}
              </div>
              
              <Button 
                type="submit"
                disabled={isGenerating}
                className="w-full sm:w-auto bg-blue-600 hover:bg-blue-700 text-white flex items-center justify-center gap-2 p-3 h-auto rounded-lg text-base font-medium"
              >
                <Beaker size={20} />
                {isGenerating ? 'Generating Molecule...' : 'Generate Molecule'}
              </Button>
            </form>
          </div>
          
          <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl">
            <h2 className="text-white text-2xl font-bold mb-4 flex items-center gap-2">
              <Atom size={24} className="text-blue-400" />
              Generated Molecules
            </h2>
            <p className="text-gray-300">No molecules generated yet. Use the form above to generate new molecules.</p>
          </div>
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default MoleculeGenerator;
