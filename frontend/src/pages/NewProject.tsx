
import React from 'react';
import { useNavigate } from 'react-router-dom';
import { useForm } from 'react-hook-form';
import { useToast } from '@/hooks/use-toast';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Textarea } from '@/components/ui/textarea';
import { 
  Form, 
  FormControl, 
  FormField, 
  FormItem, 
  FormLabel, 
  FormMessage 
} from '@/components/ui/form';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { Checkbox } from '@/components/ui/checkbox';
import { Label } from '@/components/ui/label';
import { ArrowLeft } from 'lucide-react';

// Mock existing project names for validation
const existingProjectNames = ["Project Alpha", "Project Beta", "Project Gamma", "Project Delta", "Project Epsilon"];

// Form schema type
type NewProjectFormValues = {
  name: string;
  description: string;
  targetProtein: string;
  status: string;
  therapeuticArea: string;
  collaborators: string[];
  properties: {
    solubility: string;
    bindingAffinity: string;
    toxicity: string;
    bioavailability: string;
  }
};

const NewProject: React.FC = () => {
  const navigate = useNavigate();
  const { toast } = useToast();
  const [isSubmitting, setIsSubmitting] = React.useState(false);
  
  // Initialize form with react-hook-form
  const form = useForm<NewProjectFormValues>({
    defaultValues: {
      name: '',
      description: '',
      targetProtein: 'None',
      status: 'Planning',
      therapeuticArea: 'None',
      collaborators: [],
      properties: {
        solubility: 'Medium',
        bindingAffinity: 'Medium',
        toxicity: 'Low',
        bioavailability: 'Medium',
      }
    }
  });

  // Selected collaborators state
  const [selectedCollaborators, setSelectedCollaborators] = React.useState<string[]>([]);

  // Handle form submission
  const onSubmit = (data: NewProjectFormValues) => {
    // Check if project name already exists
    if (existingProjectNames.includes(data.name)) {
      toast({
        title: "Error",
        description: "Project name already exists. Please choose a different name.",
        variant: "destructive"
      });
      return;
    }

    // Validate description length
    if (data.description.length > 500) {
      toast({
        title: "Error",
        description: "Description cannot exceed 500 characters.",
        variant: "destructive"
      });
      return;
    }

    setIsSubmitting(true);

    // In a real app, you would save the project to the database here
    // For now, we'll simulate a save with a timeout
    setTimeout(() => {
      // Create a new project object that would be saved
      const newProject = {
        id: existingProjectNames.length + 1,
        name: data.name,
        description: data.description,
        targetProtein: data.targetProtein,
        status: data.status,
        therapeuticArea: data.therapeuticArea,
        collaborators: selectedCollaborators,
        desiredProperties: data.properties,
        molecules: 0,
        lastUpdated: new Date().toISOString().split('T')[0]
      };
      
      console.log('New project created:', newProject);
      
      setIsSubmitting(false);
      toast({
        title: "Success",
        description: "Project created successfully."
      });
      navigate('/projects');
    }, 1000);
  };

  // Handle collaborator selection
  const toggleCollaborator = (collaborator: string) => {
    setSelectedCollaborators(prev => 
      prev.includes(collaborator)
        ? prev.filter(c => c !== collaborator)
        : [...prev, collaborator]
    );
  };

  return (
    <div className="flex flex-col min-h-screen bg-gray-900">
      <Navbar />
      
      <div className="flex flex-1">
        <Sidebar />
        
        <main className="flex-1 p-6">
          <div className="mb-6">
            <Button
              variant="outline"
              className="text-gray-300 border-gray-700 hover:bg-gray-800"
              onClick={() => navigate('/projects')}
            >
              <ArrowLeft size={16} className="mr-2" />
              Back to Projects
            </Button>
          </div>

          <h1 className="text-white text-3xl font-bold mb-8">Create New Project</h1>
          
          <div className="bg-gray-800/50 border border-gray-700 rounded-xl p-6 max-w-3xl">
            <Form {...form}>
              <form onSubmit={form.handleSubmit(onSubmit)} className="space-y-6">
                <FormField
                  control={form.control}
                  name="name"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel className="text-white">Project Name *</FormLabel>
                      <FormControl>
                        <Input 
                          placeholder="e.g., Project Zeta" 
                          className="bg-gray-700/50 border-gray-600 text-white"
                          {...field} 
                        />
                      </FormControl>
                      <FormMessage />
                    </FormItem>
                  )}
                />
                
                <FormField
                  control={form.control}
                  name="description"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel className="text-white">Description</FormLabel>
                      <FormControl>
                        <Textarea 
                          placeholder="Describe the goals and scope of this project..." 
                          className="min-h-[100px] bg-gray-700/50 border-gray-600 text-white"
                          {...field} 
                        />
                      </FormControl>
                      <p className="text-xs text-gray-400 mt-1">
                        {field.value.length}/500 characters
                      </p>
                      <FormMessage />
                    </FormItem>
                  )}
                />
                
                <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                  <FormField
                    control={form.control}
                    name="targetProtein"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-white">Target Protein</FormLabel>
                        <Select 
                          onValueChange={field.onChange} 
                          defaultValue={field.value}
                        >
                          <FormControl>
                            <SelectTrigger className="bg-gray-700/50 border-gray-600 text-white">
                              <SelectValue placeholder="Select a target protein" />
                            </SelectTrigger>
                          </FormControl>
                          <SelectContent className="bg-gray-800 border-gray-700">
                            <SelectItem value="None">None</SelectItem>
                            <SelectItem value="EGFR">EGFR</SelectItem>
                            <SelectItem value="HER2">HER2</SelectItem>
                            <SelectItem value="Dopamine Receptors">Dopamine Receptors</SelectItem>
                            <SelectItem value="Serotonin Receptors">Serotonin Receptors</SelectItem>
                            <SelectItem value="Histamine Receptors">Histamine Receptors</SelectItem>
                          </SelectContent>
                        </Select>
                        <FormMessage />
                      </FormItem>
                    )}
                  />
                  
                  <FormField
                    control={form.control}
                    name="status"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-white">Status *</FormLabel>
                        <Select 
                          onValueChange={field.onChange} 
                          defaultValue={field.value}
                        >
                          <FormControl>
                            <SelectTrigger className="bg-gray-700/50 border-gray-600 text-white">
                              <SelectValue placeholder="Select a status" />
                            </SelectTrigger>
                          </FormControl>
                          <SelectContent className="bg-gray-800 border-gray-700">
                            <SelectItem value="Planning">Planning</SelectItem>
                            <SelectItem value="In Progress">In Progress</SelectItem>
                            <SelectItem value="On Hold">On Hold</SelectItem>
                            <SelectItem value="Completed">Completed</SelectItem>
                          </SelectContent>
                        </Select>
                        <FormMessage />
                      </FormItem>
                    )}
                  />
                  
                  <FormField
                    control={form.control}
                    name="therapeuticArea"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-white">Therapeutic Area</FormLabel>
                        <Select 
                          onValueChange={field.onChange} 
                          defaultValue={field.value}
                        >
                          <FormControl>
                            <SelectTrigger className="bg-gray-700/50 border-gray-600 text-white">
                              <SelectValue placeholder="Select a therapeutic area" />
                            </SelectTrigger>
                          </FormControl>
                          <SelectContent className="bg-gray-800 border-gray-700">
                            <SelectItem value="None">None</SelectItem>
                            <SelectItem value="Oncology">Oncology</SelectItem>
                            <SelectItem value="Neurology">Neurology</SelectItem>
                            <SelectItem value="Cardiology">Cardiology</SelectItem>
                            <SelectItem value="Infectious Diseases">Infectious Diseases</SelectItem>
                            <SelectItem value="Metabolic Disorders">Metabolic Disorders</SelectItem>
                          </SelectContent>
                        </Select>
                        <FormMessage />
                      </FormItem>
                    )}
                  />
                </div>
                
                <div>
                  <Label className="text-white mb-2 block">Collaborators</Label>
                  <div className="space-y-2">
                    {["Dr. John Smith", "Dr. Jane Doe", "Dr. Emily Brown"].map((collaborator) => (
                      <div key={collaborator} className="flex items-center space-x-2">
                        <Checkbox 
                          id={`collaborator-${collaborator}`}
                          checked={selectedCollaborators.includes(collaborator)}
                          onCheckedChange={() => toggleCollaborator(collaborator)}
                        />
                        <label
                          htmlFor={`collaborator-${collaborator}`}
                          className="text-sm font-medium text-gray-300 leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                        >
                          {collaborator}
                        </label>
                      </div>
                    ))}
                  </div>
                </div>
                
                <div>
                  <Label className="text-white mb-3 block">Desired Molecule Properties</Label>
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div>
                      <Label className="text-gray-400 text-sm mb-1 block">Solubility</Label>
                      <Select 
                        onValueChange={(val) => form.setValue('properties.solubility', val)} 
                        defaultValue={form.getValues('properties.solubility')}
                      >
                        <SelectTrigger className="bg-gray-700/50 border-gray-600 text-white">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent className="bg-gray-800 border-gray-700">
                          <SelectItem value="High">High</SelectItem>
                          <SelectItem value="Medium">Medium</SelectItem>
                          <SelectItem value="Low">Low</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    
                    <div>
                      <Label className="text-gray-400 text-sm mb-1 block">Binding Affinity</Label>
                      <Select 
                        onValueChange={(val) => form.setValue('properties.bindingAffinity', val)} 
                        defaultValue={form.getValues('properties.bindingAffinity')}
                      >
                        <SelectTrigger className="bg-gray-700/50 border-gray-600 text-white">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent className="bg-gray-800 border-gray-700">
                          <SelectItem value="High">High</SelectItem>
                          <SelectItem value="Medium">Medium</SelectItem>
                          <SelectItem value="Low">Low</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    
                    <div>
                      <Label className="text-gray-400 text-sm mb-1 block">Toxicity</Label>
                      <Select 
                        onValueChange={(val) => form.setValue('properties.toxicity', val)} 
                        defaultValue={form.getValues('properties.toxicity')}
                      >
                        <SelectTrigger className="bg-gray-700/50 border-gray-600 text-white">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent className="bg-gray-800 border-gray-700">
                          <SelectItem value="Low">Low</SelectItem>
                          <SelectItem value="Medium">Medium</SelectItem>
                          <SelectItem value="High">High</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    
                    <div>
                      <Label className="text-gray-400 text-sm mb-1 block">Bioavailability</Label>
                      <Select 
                        onValueChange={(val) => form.setValue('properties.bioavailability', val)} 
                        defaultValue={form.getValues('properties.bioavailability')}
                      >
                        <SelectTrigger className="bg-gray-700/50 border-gray-600 text-white">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent className="bg-gray-800 border-gray-700">
                          <SelectItem value="High">High</SelectItem>
                          <SelectItem value="Medium">Medium</SelectItem>
                          <SelectItem value="Low">Low</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                  </div>
                </div>
                
                <div className="flex gap-4 pt-2">
                  <Button 
                    type="submit" 
                    className="bg-blue-500 hover:bg-blue-600 text-white" 
                    disabled={isSubmitting}
                  >
                    {isSubmitting ? "Creating..." : "Create Project"}
                  </Button>
                  <Button 
                    type="button" 
                    variant="outline" 
                    className="text-gray-300 border-gray-700 hover:bg-gray-800"
                    onClick={() => navigate('/projects')}
                  >
                    Cancel
                  </Button>
                </div>
              </form>
            </Form>
          </div>
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default NewProject;
